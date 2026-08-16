# Modal setup for Proto tools

Why this exists: `proto-tools` cannot dispatch tools locally on Windows — it raises
`Unsupported operating system: Windows` from
`proto_tools/utils/tool_instance.py::_ensure_micromamba`. Setting `device="modal"` runs the
tool in a remote Linux container instead, which is the only working dispatch path on this
machine.

Env used throughout: conda env `eab-aptamer`
(`C:\Users\christopher.brenden\AppData\Local\anaconda3\envs\eab-aptamer`). `modal` 1.5.4 is
already installed there as a proto-tools dependency — do not `pip install modal` again.

## Current state (checked 2026-08-15)

- `modal` 1.5.4 present in `eab-aptamer`.
- **Not authenticated.** No `~/.modal.toml`, no `MODAL_TOKEN_ID` / `MODAL_TOKEN_SECRET`,
  `modal profile list` is empty, and `modal token info` returns "Token missing."
- Step 1 below is interactive (opens a browser) and must be run by you.

## Credits

Modal gives new accounts $30 of free trial credit automatically. The **hackathon** credits
are separate: the re:Agent organizers hand out a redemption link at the event (check the
event Slack/Discord, the welcome email, or ask at the Modal booth). Redeem it in the Modal
dashboard under Settings -> Billing after creating your account. Neither the link nor the
credits are needed for the smoke test below — ViennaRNA is a CPU tool and costs cents.

## Setup, in order

```bash
conda activate eab-aptamer

# 1. Authenticate. INTERACTIVE — opens a browser. Run this yourself.
modal setup                      # writes ~/.modal.toml
modal token info                 # verify: prints token id + workspace

# 2. Create the environment proto-tools defaults to.
modal environment create proto-env

# 3. See what can be deployed, then deploy only what you need.
proto-tools deploy --list
proto-tools deploy --apps viennarna --env proto-env

# 4. Confirm it is live.
proto-tools deploy --status --apps viennarna --env proto-env
```

Notes:

- Each tool is its own Modal app; deploy them one at a time. `--apps all` deploys
  everything and costs real money.
- Deploying and caching model weights on a Modal volume costs money (a one-time build cost
  per tool, plus minimal ongoing storage). Remove unused weights/deployments from the Modal
  dashboard.
- `--apps <slug> --test` deploys and then smoke-tests in one shot.
- Deploys are occasionally flaky on third-party download links; retry once before
  investigating.
- Optional: `export PROTO_MODAL_SCALEDOWN_WINDOW=300` keeps containers warm for 5 min
  (faster repeat calls, more idle cost). Default is 30 s.
- `MODAL_ENVIRONMENT=proto-env` can replace `--env proto-env` everywhere.

## Dispatching to Modal

Set `device="modal"` on the tool's **config** object. Nothing else changes.

```python
from proto_tools import run_viennarna, ViennaRNAInput, ViennaRNAConfig

out = run_viennarna(
    ViennaRNAInput(sequences=[seq]),
    ViennaRNAConfig(device="modal"),
)
```

The tool must be deployed (step 3) before `device="modal"` will work.

## Smoke test

Verifies the Modal dispatch path and cross-validates the project's thermodynamics against
an independent execution of ViennaRNA. Run after step 4.

```bash
python -c "
from proto_tools import run_viennarna, ViennaRNAInput, ViennaRNAConfig
seq = 'TGTGTATAAGTCACGAGGTTTGAATAAGAACCATCGGCGCCAACAAAACATTCAAACAGAAATCTACTAGTCAC'
r = run_viennarna(
    ViennaRNAInput(sequences=[seq]),
    ViennaRNAConfig(device='modal', use_dna_params=True, temperature=37.0, no_lonely_pairs=True),
).results[0]
print(r.structure)
print(f'{r.mfe:.2f} kcal/mol')
"
```

Expected output — the locally-computed reference from the standalone ViennaRNA package for
IL-6-9805 (the current parent aptamer):

```
.(((.......)))...((((((((........................)))))))).................
-4.10 kcal/mol
```

An exact match means the Modal path works and the two engines agree. Any difference in the
structure string or MFE is a real discrepancy — report it rather than rounding it away.
