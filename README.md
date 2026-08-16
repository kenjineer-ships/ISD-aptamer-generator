# Structure-switching IL-6 aptamers for electrochemical aptamer-based sensors (E-AB)

Dry-lab pipeline that produces a **ranked shortlist of DNA aptamer sequences** predicted to
both bind IL-6 and structure-switch well enough to transduce signal in an E-AB biosensor.

re:Agent Hackathon, 15–16 August 2026.

## Setup

**Windows users — do this first.** The `proto-language` repo contains example files with very
long names. pip clones it into a deep temp directory, and the combined path exceeds Windows'
260-character `MAX_PATH` limit, so `git checkout` fails mid-clone and the whole env build
aborts with `CondaEnvException: Pip failed`. One-time, per-user, no admin needed:

```bash
git config --global core.longpaths true
```

Then, on any platform:

```bash
git clone <this repo>
cd "ReAgent Hackathon"
conda env create -f environment.yml
conda activate eab-aptamer
```

Then check it works:

```bash
python aptamer/build_parents.py
```

You should see four IL-6 aptamers with their secondary structures, MFE values, and a
fractional-occupancy table.

### Platform note

The Windows limitation is narrower than "Proto doesn't work here". It applies only to
**tool dispatch**, not to the language layer.

`proto-tools` raises `Unsupported operating system: Windows` from
`proto_tools/utils/tool_instance.py::_ensure_micromamba`, which is reached only when a tool is
actually dispatched. Importing `proto_tools` and `proto_language` succeeds, and the language
layer — `Segment`, `Construct`, `Constraint`, generators, and the MCMC/GA/rejection-sampling
optimizers — is pure Python and dispatches nothing.

Because a custom constraint is just a function returning `ConstraintOutput(score=...)`, a
constraint that calls the standalone `RNA` and `primer3` packages never touches the micromamba
layer. **The entire thermodynamic switch-design loop therefore runs natively on Windows.**
Verified by `aptamer/test_proto_language.py`, which drives a proto-language MCMC over a
linker + displacement-strand construct with a ViennaRNA-backed constraint.

| | Windows | Linux / macOS |
|---|---|---|
| ViennaRNA / primer3 thermodynamics | native | native |
| proto-language optimizers + custom constraints | **native** | native |
| Built-in tool-backed constraints (`boltz2-binding-strength`, `af3-*`, `gyration-radius`, Proto's own `viennarna` wrapper) | `device="modal"` or WSL2 | native |
| GPU co-folding (Boltz2, AF3, OpenDDE, Protenix) | `device="modal"` or WSL2 | native |

### proto-language gotcha

`Generator.assign()` **replaces** the generator's segment list rather than appending, and
passing several segments at once *ties* them (shared generated values, identical length
required). Use **one generator per independently-varying segment**:

```python
gen_linker = RandomNucleotideGenerator(cfg); gen_linker.assign(linker)
gen_ds     = RandomNucleotideGenerator(cfg); gen_ds.assign(ds)
MCMCOptimizer(..., generators=[gen_linker, gen_ds], ...)
```

Calling `assign()` twice on one generator silently orphans the first segment, and the
optimizer then fails with `references segment 'linker' which has no populated sequence and no
generator assigned`.

## Approach

The pipeline deliberately **does not** rank candidates by comparing a predicted unbound 3D
structure to a predicted bound 3D structure. That approach is not supported by the current
state of the field:

- CASP16 target D1273 was a 27-nt DNA aptamer. Across 107 models from 22 groups, none
  achieved RMSD below 10 Å (range 10.16–27.79 Å); ~70% of predicted base pairs were
  non-canonical.
- AlphaFold3 on aptamers averages 1.45 Å RMSD for structures deposited *before* its training
  cutoff versus 6.40 Å for those after — novel sequences fall in the second regime.
- Boltz-2's own paper notes the models "often fail to capture large conformational changes,
  such as those that can be induced by binding" — exactly the quantity of interest.

Instead we engineer an explicit **intramolecular strand displacement (ISD)** switch:

```
5'-[ aptamer ]--[ poly-T linker ]--[ displacement strand ]-3'
```

The displacement strand (DS) hybridises to the aptamer's binding region; target binding
displaces it. Switching becomes a designed three-state equilibrium whose energetics are
computable with nearest-neighbour thermodynamics, rather than an emergent property we hope a
diffusion model resolves.

```
State A: DS hybridised (closed)      State B: DS open      State C: target bound

dG_switch  = G(A) - G(B)
K_switch   = exp(-dG_switch / RT)
KD_app     = KD_parent * (1 + K_switch)     <- switching costs affinity
gain_proxy = f_closed(no target) - f_closed(+target)
```

Signal gain in an E-AB is the *difference in electron-transfer rates* between unbound and
bound states — not a 5'-to-3' distance — so `gain_proxy` (a population shift) is the
appropriate objective.

## Layout

```
environment.yml           conda environment
dashboard.html            output: self-contained browser dashboard (open by double-click)
docs/modal-setup.md       Modal remote-compute setup (needed for the GPU stages on Windows)
aptamer/
  build_parents.py        reconstruct + fold the four published IL-6 aptamers
  parents.json            output: sequences, structures, MFE, SPR kinetics
  1ALU.pdb                human IL-6 crystal structure, 1.9 Å (receptor input)
  switch_library.py       enumerate + rank ISD switch constructs -> switches.csv
  mismatch_tune.py        refine short-tether constructs with single DS mismatches
  mismatches.csv          output: passing single-mismatch variants (dashboard view 9)
  negative_controls.py    do composition-matched controls separate from designed DS?
  negative_controls.csv   output: designed + scrambled/reversed/foreign arms (view 10)
  cofold/                 co-folded aptamer-IL-6 complexes; dashboard views 11-12 and the
                          3D panel read it (read-only) -- see cofold/README.md
  benchmark.py            validate the model against published ISD measurements
  test_proto_language.py  smoke test: proto-language MCMC + ViennaRNA constraint on Windows
  dashboard.py            read the outputs above -> dashboard.html (all data inlined)
```

Verify your setup, then run the design:

```bash
python aptamer/build_parents.py         # four aptamers, structures, occupancy table
python aptamer/test_proto_language.py   # ends with "PROTO-LANGUAGE ON WINDOWS: OK"
python aptamer/switch_library.py        # ranked shortlist -> aptamer/switches.csv
python aptamer/dashboard.py             # regenerate ./dashboard.html, then open it
```

`dashboard.html` inlines every value it shows, so it works from `file://` with no server and
no build step — never hand-edit it, regenerate it. Missing pipeline outputs render as
"not yet generated" rather than failing. Two views prefer a CDN script: the 3D structure
panel loads 3Dmol.js and genuinely needs internet, and the secondary-structure panel loads
fornac for the real 2D layout but falls back to the built-in base-pair arc diagram (and says
so in the panel) when fornac is unreachable. The tables, occupancy curve, scatter, heatmap,
coverage track and the two co-fold views need no network at all. Only five of the 32 predicted
complexes are inlined into the 3D panel (about 1.1 MB of the page); the rest stay on disk in
`aptamer/cofold/structures/`.

### Switch design target

The DS-stabilisation energy is **derived from the affinity budget, not chosen**:

```
dG_switch = G(closed) - G(open)          G(open) = MFE(aptamer) + MFE(linker+DS)
K_closed  = exp(-dG_switch / RT)         KD_app  = KD_parent * (1 + K_closed)

filters:  KD_app <= 1 uM  (affinity budget)
          closed fraction >= 75%  (the switch must actually close)
```

Filter on those two quantities directly, **not** on a dG window. The acceptable dG depends on
the parent's K_D, so a hardcoded window leaks over-budget constructs the moment the parent
changes — which it did when the parent moved from IL-6-7326 to IL-6-9805.

**Parent choice.** The parent must be one whose published K_D was measured on the *full-length*
sequence. IL-6-7326 is tighter (27 nM) but that number belongs to the truncated 7326.1, so
pairing it with the full 74-mer optimises against a different molecule. `switch_library.py`
asserts this.

### Benchmark: does the model reproduce measured reality?

`python aptamer/benchmark.py`

Ground truth is Wilson, Hariri, Thompson, Eisenstein & Soh, *Nat Commun* (2019), PMC6838323,
Supplementary Table 2: the closing equilibrium constant **K_Q** measured on an **ATP** aptamer
across DS and loop lengths, plus their fitted design law `d(log10 K_Q)/d(L_DS) = 0.826 ± 0.157`.
K_Q is the same quantity this pipeline calls `K_closed`.

Run either parent: `python aptamer/benchmark.py --parent IL-6-9805`

| Test | IL-6-7326.1 (45 nt) | IL-6-9805 (74 nt) |
|---|---|---|
| Slope of log K vs DS length (published 0.826 ± 0.157) | +0.516 **out of range** | +0.627 **out of range** |
| Loop-length dependence | correct sign at every DS, within 15–25% of measured | " |
| Spearman vs 12 measured grid points | +0.972 | **+1.000** |

**Honest reading: the model ranks correctly but under-predicts DS-length sensitivity by
25–40%.** Rank agreement is near-perfect and the loop-length knob is quantitatively close, but
adding a DS base changes K less in this model than Wilson measured. Practical consequence: the
K_D,app spread across DS lengths is *compressed* — the real affinity penalty for a long DS is
worse than the table says, so prefer short DS more strongly than the numbers imply.

Two earlier claims in this file were wrong and are corrected above:

- An earlier version reported **+0.818 on the 74-mer, "in range"**. That used
  `K_closed = exp(-dG/RT)`, which overstates K_closed by exactly 1 (see `k_closed_from`). The
  apparent pass was an artifact of that bias.
- It also reported the validation as strongly **parent-dependent** (+0.818 vs +0.350). With the
  unbiased estimator the gap shrinks to +0.627 vs +0.516. Most of the apparent
  parent-dependence was the biased estimator interacting with weak switches, not the parent.

A `--engaged-only` aggregation was also tried and made things worse (slope −0.166, Spearman
−0.448): conditioning on windows where the DS wins is selection bias, since more windows
qualify as DS length grows, so the median's composition shifts with the variable being
regressed on. The flag is kept for inspection; the default is all windows.

### Negative controls

`python aptamer/negative_controls.py`

Appending *any* sequence lowers a construct's free energy, so a negative dG_switch proves
nothing on its own. Each designed DS is compared against three controls matched on length and
base composition: `scrambled` (shuffled), `reversed` (reversed, not complemented), and
`foreign` (reverse complement of a window from a *shuffled* aptamer — a genuine duplex-former
that simply isn't complementary to ours; the strictest of the three).

| arm | median dG | median closed | median engagement | pass rate |
|---|---|---|---|---|
| designed | −1.50 | 91.3% | **0.83** | 100%* |
| scrambled | −0.25 | 33.4% | 0.02 | 1% |
| reversed | −0.26 | 33.7% | 0.03 | 4% |
| foreign | −0.30 | 38.9% | 0.05 | 0% |

Zero of 240 controls reach the designed median engagement. The score is specific to
complementarity with this aptamer, not to DNA-likeness.

\* The designed arm's 100% is partly circular — those constructs were drawn from
`switches.csv`, which already passed these filters. The non-circular evidence is the
engagement gap and the *control* pass rates, since controls were never pre-selected.

Note the controls still show dG_switch of −0.25 to −0.30 kcal/mol. That residual is the
trivial signal that raw dG_switch alone would have admitted; the engagement filter is what
excludes it. Treat 1–4% as the false-positive floor.

What changed the numbers, in order of effect:

| Fix | Slope (45 nt) | Loop test | Spearman |
|---|---|---|---|
| MFE, `K = exp(-dG/RT)` | +0.350 | wrong sign at DS 6–7, 2–3× too small | +0.846 |
| + ensemble dG instead of MFE | +0.434 | — | +0.972 |
| + `K = exp(-dG/RT) - 1` | +0.516 | correct everywhere, within 15–25% | +0.972 |

**Consequence for the design.** There is a real tension to resolve:

| Parent | Construct length | Model validated? | K_D |
|---|---|---|---|
| IL-6-9805 (74 nt, full) | 83–128 nt — too long for a good E-AB | yes | 40.5 nM |
| IL-6-7326.1 (45 nt, truncated) | 54–100 nt — buildable | **no** | 27 nM (tightest) |

The pipeline currently uses the truncated parent, i.e. the better device with the weaker
validation. Switch `PARENT` in `switch_library.py` to trade back.

**Read the caveats before quoting these numbers.**

- The Spearman of +0.993 is weaker evidence than it looks. Both predicted and measured K rise
  steeply with DS length, so any model with the right DS trend scores high on a 12-point grid.
- **Loop-length sensitivity is under-predicted 2–3×** (−0.027 vs −0.078 per nt at DS 9). The
  nearest-neighbour model treats a poly-T linker as near-energy-neutral and its loop-entropy
  penalty is probably under-parameterised at these lengths. Consequence: the linker knob is
  likely *more* powerful in the wet lab than `switches.csv` implies. Trust the ordering across
  linker lengths, not the absolute spacing.
- This validates the **design laws, not per-construct prediction**. Their DS sequences are not
  published in reachable form, so K_Q could not be recomputed for their exact constructs.
- The value `3.58` appears both in this K_Q matrix and as `k_obs` in the review's Figure 1d
  (PMC11883736), which raised the worry that the extracted table was kinetics rather than
  thermodynamics. Recomputing their slope from the matrix gives 0.82–1.03, consistent with the
  stated 0.826, which resolves it in favour of K_Q.
- The matrix was read from a summary of the supplementary spreadsheet, giving a 3×5 sub-block
  (loops 23/25/33 × DS 9/8/7/6). PMC serves a bot-mitigation interstitial instead of the
  `.xlsx`, so the full grid was unavailable. DS=5 was excluded as unfittable noise
  (21.4 / 0.0379 / 9.2, non-monotonic in both directions).

Two traps worth not re-falling into:

- A plausible-sounding target like −7.5 kcal/mol is catastrophic: `K_closed ~ 2e5`, so
  `KD_app ~ 5 mM` and the switch never opens.
- `MFE(construct) - MFE(aptamer)` is the **wrong** estimator — it absorbs the DS's own
  folding. The open state must include `MFE(linker+DS)`.

Ranking cannot be by `KD_app`, which is monotone in `dG_switch` and so just returns the window
edge; a (K_D, gain) Pareto front is degenerate for the same reason. Constructs are instead
discriminated by how much **randomised** library sequence the DS covers — selection only
varied those positions, so the binding site almost certainly involves them — then by
engagement, then by short DS length to limit cross-hybridisation.

## Parent aptamers

Reconstructed from the Neomer library template plus per-candidate random-module identities.
SPR kinetics as published; none showed measurable binding to human serum albumin.

| Aptamer | K_D | k_off | Residence time | MFE @37 °C |
|---|---|---|---|---|
| IL-6-7326 | 27.0 nM | 1.50e-3 /s | 667 s | −5.70 |
| IL-6-6449 | 32.1 nM | 1.52e-3 /s | 658 s | −5.90 |
| IL-6-9805 | 40.5 nM | 1.77e-3 /s | 565 s | −4.10 |
| IL-6-4202 | 73.2 nM | 2.98e-3 /s | 336 s | −3.40 |

**Verify before ordering oligos:**

1. `IL-6-6449`'s third random block reads as `GG` (2 nt) where the template expects 3 nt —
   one nucleotide is uncertain. Prefer 7326 or 9805.
2. The source text states a 73-nt library, but the printed template parses to 74. A 1-nt
   ambiguity in the fixed regions changes the predicted fold. Check the supplementary
   material or the IDT order sheet.

## Known constraint: the affinity gap

Wound fluid is a far better matrix than serum for IL-6 (100–1000× higher), but the parent
aptamers are still weak relative to physiological concentrations:

```
IL-6-7326 (KD 27 nM) fractional occupancy in wound fluid
  post-surgical WHD (median)     987 pg/mL =   42 pM ->  0.2% bound
  inflammatory (median)         4964 pg/mL =  209 pM ->  0.8% bound
  infected (median)             5883 pg/mL =  248 pM ->  0.9% bound
  max observed                135500 pg/mL = 5713 pM -> 17.5% bound
```

Adding a displacement strand raises K_D further.

**Decision (noted, proceeding):** the numbers above are bulk-equivalent concentrations from
diluted samples. The skin model returned up to **5 ng/mL from a likely heavily diluted
sample**, so true undiluted levels are higher, and a high-spatiotemporal-resolution electrode
array samples *local* concentration at the sensor surface rather than bulk. Both effects push
real occupancy well above the table. Revisit if measured signal-to-noise disappoints; the
fallbacks are repositioning as a threshold alarm for severe inflammation, or inserting an
affinity-maturation stage ahead of switch engineering.

This decision sets the affinity budget the switch design optimises against: `KD_app <= 1 uM`.

Residence times of 5.6–11 minutes also make these intrinsically slow for continuous
monitoring.

## Sources

- Meehan et al. *A new method for the reproducible development of aptamers (Neomers).*
  PLOS ONE (2025). doi:10.1371/journal.pone.0311497 — parent aptamers, library template,
  SPR kinetics.
- Saunders, Thompson & Soh. *Generalizable Molecular Switch Designs for In Vivo Continuous
  Biosensing.* (2025). PMC11883736 — ISD/DBS switch architectures and tuning rules.
- Bakestani et al. *Carboxylate-Terminated Electrode Surfaces Improve the Performance of
  Electrochemical Aptamer-Based Sensors.* ACS AMI (2025). doi:10.1021/acsami.4c21790 —
  E-AB gain defined by electron-transfer rate difference.
- Kretsch et al. *Functional relevance of CASP16 nucleic acid predictions.* (2025).
  PMC12412911 — limits of nucleic-acid structure prediction.
- Rembe et al. *Immunomarker profiling in human chronic wound swabs.* (2025). PMC11978031 —
  wound-fluid IL-6 concentrations.
- IL-6: UniProt P05231; structure PDB 1ALU (1.9 Å).
