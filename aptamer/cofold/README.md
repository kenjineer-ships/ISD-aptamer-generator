# Aptamer–IL-6 co-folds

Predicted structures for the IL-6 aptamer and its top-ranked switch constructs, co-folded
against human IL-6 and (for specificity) the rest of the gp130 cytokine family, using four
structure-prediction models deployed to Modal via `proto-tools`.

> **These structures are figures and a weak specificity signal, NOT ranking input.**
> The pipeline deliberately does not rank on predicted structure. See
> [Why these numbers do not rank anything](#why-these-numbers-do-not-rank-anything).

**Status: complete. 32/32 predictions succeeded**; 32 structures in `structures/`, `failures` in
`cofold_summary.json` is empty. All four models deployed and ran on H100.

## Headline findings

The headline results are negative ones, and that is the point.

**1. The four models do not agree on where the aptamer binds IL-6.**

| DNA chain | mean pairwise Jaccard | tolerant ±2 | residues all 4 models agree on | union |
|---|---|---|---|---|
| `parent45` | 0.168 | 0.236 | **0** | 34 |
| `construct1` | 0.162 | 0.324 | **0** | 82 |
| `construct2` | 0.180 | 0.302 | 1 | 69 |
| `construct3` | 0.143 | 0.211 | 1 | 65 |
| `construct4` | 0.135 | 0.217 | **0** | 70 |
| `construct5` | 0.127 | 0.251 | **0** | 74 |

Out of unions of 34–82 contact residues, the number every model agrees on is zero or one. Some
pairs share literally nothing: `boltz2` vs `opendde` has Jaccard **0.000** on `parent45`,
`construct1` and `construct4` — completely disjoint binding sites.

**2. Three of the four broadly agree; OpenDDE is the outlier.** `boltz2`, `esmfold2` and
`protenix` all place the DNA on two patches around residues **~11–30** and **~109–127**, while
`opendde` puts it on a different face (~39–53, ~150–183) and reports the *highest* confidence
while doing so. Restricted to the agreeing trio, consensus residues recur across chains:
**26/27, 30, 117, 120, 123/124**. That is a hypothesis to test experimentally, not a result.

**3. The specificity signal fails its own control.** For the bare 45-nt parent — the chain with
an actual measured Kd of 27 nM for IL-6 — Boltz-2 scores **CNTF (0.597) higher than the true
target IL-6 (0.581)**, a margin of **−0.016**. The model cannot identify the correct target among
five homologues. Every off-target number below must be read in light of that.

**4. Confidence collapses with tether length.** Adding the poly-T linker drops Boltz-2 ipTM from
0.581 (`parent45`, 45 nt) to 0.173–0.367 (constructs, 72–99 nt), and DNA-chain pTM from 0.696 to
0.36–0.49. This is largely an artefact of scoring 27–54 nt of disordered spacer, not evidence
that the constructs bind worse.

## What was run

**Phase A — cross-model agreement.** 6 DNA chains × 4 models against IL-6 = 24 predictions.

**Phase B — off-target specificity.** Boltz-2 only, against the gp130 family: IL-6 (`P05231`,
on-target), IL-11 (`P20809`), LIF (`P15018`), OSM (`P13725`), CNTF (`P26441`). Run for both
`construct1` (the top-ranked construct, as specified) and `parent45` (added as a control, because
the parent is the chain with real measured affinity and a much cleaner interface score).
8 new predictions; the two on-target runs are reused from Phase A.

| DNA chain | Length | What it is |
|---|---|---|
| `parent45` | 45 nt | the bare parent aptamer, `IL-6-7326.1`, Kd 27 nM |
| `construct1`…`construct5` | 72–99 nt | the top 5 rows of `switches.csv`, in existing pipeline rank order |

### Models and deploy outcomes

| Model | proto-tools key | Deploy | Notes |
|---|---|---|---|
| Boltz-2 | `boltz2-prediction` | 12.6 min (2nd attempt) | MSA-based (MMseqs2), diffusion decoder |
| Protenix | `protenix-prediction` | 11.1 min (3rd attempt) | open-source AlphaFold3 reimplementation, MSA-based |
| OpenDDE | `opendde-prediction` | 7.0 min | all-atom biomolecular complex prediction |
| ESMFold2 | `esmfold2-prediction` | 5.8 min | all-atom, **single-sequence (no MSA)** |

All ran on H100 in Modal environment `proto-env`, stock tool configs, `seed=42`, `device="modal"`.

### Sequences

**Protein.** Mature chains only, taken from each UniProt entry's own `Chain` feature, so the
signal peptide (and for OSM the C-terminal propeptide) is stripped. The signal peptide is cleaved
before IL-6 ever encounters an aptamer; leaving it on would have the model fold 29 hydrophobic
residues that do not exist in the mature cytokine.

| Target | UniProt | Full | Mature chain used | Length used |
|---|---|---|---|---|
| IL-6 | P05231 | 212 aa | 30–212 (signal 1–29 removed) | 183 aa |
| IL-11 | P20809 | 199 aa | 22–199 (signal 1–21 removed) | 178 aa |
| LIF | P15018 | 202 aa | 23–202 (signal 1–22 removed) | 180 aa |
| OSM | P13725 | 252 aa | 26–221 (signal 1–25 and propeptide 222–252 removed) | 196 aa |
| CNTF | P26441 | 200 aa | 1–200 (cytosolic, no signal peptide) | 200 aa |

**Residue numbering in this directory is 1-based within the mature chain**, so IL-6 residue 1 is
UniProt P05231 residue 30. Add 29 to convert to UniProt numbering.

**DNA.** Taken from the pipeline, never retyped: `parent45` from
`build_parents.parent("IL-6-7326.1")`, constructs from the `construct` column of `switches.csv`.
Because `switches.csv` is regenerated while the pipeline is being edited, it was snapshotted at
run start — `switches_snapshot.csv` and `inputs_snapshot.json` (with a SHA-256) record exactly
which version these predictions correspond to.

**Chain convention.** Chain `A` is always the protein; chain `B` is always the DNA.

## Results

### Phase A — interface confidence (ipTM) against IL-6

| DNA | nt | boltz2 | protenix | opendde | esmfold2 |
|---|---|---|---|---|---|
| `parent45` | 45 | 0.581 | 0.561 | 0.632 | 0.279 |
| `construct1` | 94 | 0.367 | 0.246 | 0.643 | 0.190 |
| `construct2` | 98 | 0.258 | 0.349 | 0.633 | 0.082 |
| `construct3` | 99 | 0.208 | 0.307 | 0.685 | 0.081 |
| `construct4` | 72 | 0.321 | 0.358 | 0.667 | 0.187 |
| `construct5` | 95 | 0.173 | 0.223 | 0.641 | 0.213 |

Note `opendde` returns a high, nearly constant ipTM (0.63–0.69) regardless of chain — including
for the chains where it disagrees with every other model about the binding site. A confidence
score that does not vary with the input is not carrying information.

### Phase A — where the DNA lands (contact patches, mature-chain numbering)

| DNA | boltz2 | esmfold2 | protenix | opendde |
|---|---|---|---|---|
| `parent45` | 15,23–27,109–124 | 14–18,23–30,109–120,127 | 11–15,26–30,36,117–124,170–178 | 39,50–53,99,150,165–170,178–181 |
| `construct1` | 1,13–16,23–33,40,109–124 | 15–30,36,49–51,74,94,112–123,137,170–178 | 23–47,53–60,65,106–110,120–127,149,154–157,167,174–181 | 20,39,49–74,90,96,118,126–130,149,154–159,167–173,178–183 |

Full per-chain patches for all six DNA chains are in
`cofold_summary.json → cross_model_agreement.per_dna.*.contact_patches_per_model`.

Outlier check (mean tolerant-Jaccard of each model against all others, averaged over chains):
`protenix 0.316`, `esmfold2 0.297`, `opendde 0.210`, `boltz2 0.204`. No model is close to
agreement with the rest; `opendde`'s disagreement is the most systematic because it is
confidently on a different face.

### Phase B — off-target margins (Boltz-2 ipTM, on-target minus off-target)

| Target | UniProt | `construct1` ipTM | margin | `parent45` ipTM | margin |
|---|---|---|---|---|---|
| **IL-6 (on-target)** | P05231 | **0.367** | — | **0.581** | — |
| IL-11 | P20809 | 0.348 | +0.019 | 0.389 | +0.192 |
| LIF | P15018 | 0.086 | +0.281 | 0.397 | +0.184 |
| OSM | P13725 | 0.236 | +0.131 | 0.381 | +0.201 |
| CNTF | P26441 | 0.222 | +0.145 | **0.597** | **−0.016** |

Read this as a negative result, not a specificity claim:

- For `construct1`, IL-11 is separated from IL-6 by **+0.019 ipTM** — noise. The model does not
  distinguish the on-target from the closest family member.
- For `parent45`, **CNTF outscores IL-6**. The one chain here with a measured 27 nM Kd for IL-6 is
  assigned its *highest* interface confidence against the wrong cytokine.
- The two panels do not even agree with each other on the ordering of off-targets (LIF is the
  best-discriminated target for `construct1` and the worst for `parent45`).

The honest summary: **Boltz-2 ipTM does not resolve gp130-family specificity for this aptamer.**
That is useful to know — it rules the metric out as a selection filter — but it is not evidence
that the aptamer is or is not selective.

## Files

| Path | Contents |
|---|---|
| `structures/<model>__<dna>__<target>.cif` | predicted complex, mmCIF (32 files) |
| `results/<model>__<dna>__<target>.json` | one prediction: every metric, plus its contact epitope |
| `cofold_summary.json` | all predictions + cross-model agreement + off-target margins |
| `inputs_snapshot.json`, `switches_snapshot.csv` | exact inputs, with SHA-256 |
| `targets.json` | mature target sequences as fetched from UniProt |
| `run_cofold.py` | driver; resumable, one job at a time |
| `summarize.py` | builds `cofold_summary.json` from `results/` |

Reproduce with:

```bash
python run_cofold.py --phase A --models boltz2,protenix,opendde,esmfold2
python run_cofold.py --phase B --models boltz2 --best-construct construct1
python run_cofold.py --phase B --models boltz2 --best-construct parent45
python summarize.py
```

`device="modal"` on the tool config is what routes a call to the deployed app. Without it
proto-tools runs the model locally, which on Windows fails with
`Unsupported operating system: Windows (arch: AMD64)`. Deploying is not sufficient.

## What the metrics mean

Recorded verbatim per prediction under `metrics`; not every model emits every field.

| Metric | Range | Better | Meaning |
|---|---|---|---|
| `iptm` | 0–1 | higher | **interface** predicted TM-score — confidence in the *relative placement* of the two chains. The number that matters for "does the aptamer bind here", and the one to distrust most. |
| `ptm` | 0–1 | higher | predicted TM-score for the whole complex, dominated by the protein's own fold |
| `chains_ptm` | 0–1 | higher | per-chain pTM; for chain `B` this is confidence in the DNA's *internal* fold |
| `pair_chains_iptm` | 0–1 | higher | ipTM per chain pair; the `A`/`B` off-diagonal entry is the protein–DNA interface |
| `complex_plddt` | 0–1 | higher | mean per-atom local confidence over the complex |
| `complex_iplddt` | 0–1 | higher | pLDDT restricted to interface atoms |
| `complex_pde` / `complex_ipde` | Å | lower | predicted distance error, all atoms / interface atoms |
| `avg_pae` / `pae` | 0–32 Å | lower | predicted aligned error; cross-chain blocks are the interface uncertainty |
| `confidence_score` | 0–1 | higher | each model's own blended headline score. **Not comparable across models.** |

**`epitope`** in each record is not a model output — it is measured from the returned coordinates:
every IL-6 residue with at least one heavy atom within **4.0 Å** of any DNA atom
(`Structure.interface_contact_residues`). This is the quantity compared across models, because it
is a geometric fact about the structure rather than a self-reported score.

**`cross_model_agreement`** reports, per DNA chain, the pairwise Jaccard overlap of those epitope
residue sets, both residue-exact and with ±2 residues of tolerance (`..._tolerant_pm2`; two models
can pick the same surface patch and still score near zero if their contact lists are offset by a
residue or two). High overlap means independent models put the DNA in the same place. That is the
only confidence statement here worth anything — and as shown above, it is low.

## Why these numbers do not rank anything

Predicted protein–DNA-aptamer structures are not currently accurate enough to rank binders, and
the published failure modes are specific and severe:

- **CASP16** included a 27-nt DNA aptamer target. **0 of 107** submitted models achieved under
  10 Å RMSD. Not "poor accuracy" — no group solved it.
- **AlphaFold3 memorises aptamers.** It averages **1.45 Å** on aptamer complexes published before
  its training cutoff and **6.40 Å** on ones published after. A good-looking prediction on a known
  aptamer is a lookup, not a generalisation, and every construct here is novel.
- **Interfaces are mostly wrong even when the fold is right.** Only **35.6%** of AF3
  protein–nucleic-acid predictions recover more than half the native contacts.

Consequences for reading this directory:

1. **ipTM is not an affinity and not a rank.** Do not order constructs by it. The pipeline's
   ranking comes from `switches.csv` (thermodynamics and measured Kd), and these predictions were
   run *after* that ranking, not into it.
2. **Cross-model agreement is the signal — and here it is weak.** Four models disagreeing about
   where the DNA sits means the binding mode is unresolved, regardless of how confident any one of
   them is. Even the three that broadly converge share training data and biases, so their agreement
   is correlated evidence, not independent replication.
3. **The off-target margins are the most defensible output in principle** — same model, same DNA,
   same protocol, five homologous cytokines, so systematic errors partly cancel. In practice the
   `parent45` control (CNTF > IL-6) shows the metric failing on the one case where the right answer
   is known. Do not quote these margins as specificity.
4. **The long poly-T linkers are not real structure.** Constructs 1–5 carry 27–54 nt of poly-T
   tether. Single-stranded poly-T is a flexible, largely unstructured spacer with no defined
   conformation, so whatever the models draw there is arbitrary. It also drags down whole-complex
   scores purely by adding disordered residues — which is most of the `parent45`-vs-construct ipTM
   gap in the table above. Compare across constructs, never parent-vs-construct.
5. **Single-sequence vs MSA models are not equivalent.** ESMFold2 sees no MSA; the other three do.
   Some of the between-model spread is that difference, not genuine uncertainty about the aptamer.
   ESMFold2's uniformly low ipTM is partly this.
6. **`opendde`'s high, flat confidence is a warning sign, not a strength.** It is the most
   confident model and the one that agrees least with the others.
7. **No experimental validation is present in this directory.** Nothing here has been compared to a
   solved structure, a footprint, or a mutational scan. Use these as figures and as a hypothesis
   generator for which IL-6 surface to probe experimentally — the recurring 26/27, 30, 117, 120,
   123/124 consensus patch is the obvious first thing to mutate.

## Provenance caveats

- `switches.csv` was being regenerated by the pipeline during this run. The snapshot files pin what
  was actually used (SHA-256 in `inputs_snapshot.json`); re-running `summarize.py` against a newer
  `switches.csv` will not match. The five constructs used were rows 1–5 at snapshot time, with
  `kd_app_nM` 971.2 / 984.2 / 951.8 / 922.5 / 885.9.
- Modal initially refused to schedule H100/A100 on this workspace ("Please add a payment method to
  use H100 GPU functions"); only T4/L4/A10 were available. This was resolved mid-run. **All 32
  reported predictions ran on stock `GPU_DEFAULT` (H100) with unmodified tool definitions** — a
  temporary local patch to `gpu_profiles.py` was reverted before any of them ran, and no prediction
  in this directory came from the downgraded-GPU path.
- An earlier batch of 18 runs failed outright because dispatch went local instead of to Modal;
  those records were deleted and rerun. They are not represented here.
- One Boltz-2 and one Protenix deploy failed transiently before succeeding (a truncated CCD tarball
  download and an unexplained image-build failure respectively); both succeeded on retry with no
  configuration change.
