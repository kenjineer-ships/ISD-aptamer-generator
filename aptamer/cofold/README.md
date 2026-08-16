# Aptamer–IL-6 co-folds

Predicted structures for the IL-6 aptamer and its top-ranked switch constructs, co-folded
against human IL-6 and (for specificity) the rest of the gp130 cytokine family, using four
structure-prediction models deployed to Modal via `proto-tools`.

> **These structures are figures and a weak specificity signal, NOT ranking input.**
> The pipeline deliberately does not rank on predicted structure. See
> [Why these numbers do not rank anything](#why-these-numbers-do-not-rank-anything).

---

## STATUS: NO STRUCTURES WERE PRODUCED. ALL 18 RUNS FAILED.

Every record in `results/` has `"success": false`. There are no `.cif` or `.pdb` files in this
directory and no `cofold_summary.json`, because `summarize.py` never ran. Everything below
describes what was *attempted*, not what was obtained.

Cause: the runs executed **locally**, not on Modal, and proto-tools cannot dispatch tools on
Windows — `RuntimeError: Unsupported operating system: Windows (arch: AMD64)`, raised from
`proto_tools/utils/tool_instance.py::_ensure_micromamba`. Deploying an app to Modal is not
sufficient; the call must also be routed there with `device="modal"` on the tool's config.

Separately, `protenix` never deployed at all: all four apps request `H100/H200/A100-80GB` by
default, and Modal refuses those GPU classes until a payment method is on file — hackathon
credits alone do not unlock them. A workaround patch to `GPU_DEFAULT` (`A10:1`/`L4:1`) is
applied in site-packages, with `gpu_profiles.py.orig-backup` alongside it.

What is still useful here: `run_cofold.py`, `targets.json`, `inputs_snapshot.json` and
`switches_snapshot.csv` record exactly which sequences and targets were queued, so a rerun on
Linux, WSL2, or with `device="modal"` wired through needs no re-derivation.

---

## What was run

**Phase A — cross-model agreement.** 6 DNA chains × 4 models, each co-folded with IL-6:

| DNA chain | Length | What it is |
|---|---|---|
| `parent45` | 45 nt | the bare parent aptamer, `IL-6-7326.1` |
| `construct1`…`construct5` | 72–99 nt | the top 5 rows of `switches.csv`, in existing pipeline rank order |

**Phase B — off-target specificity.** The single best construct against the gp130 family,
`boltz2` only: IL-6 (`P05231`, on-target), IL-11 (`P20809`), LIF (`P15018`),
OSM (`P13725`), CNTF (`P26441`).

### Models

| Model | proto-tools key | Notes |
|---|---|---|
| Boltz-2 | `boltz2-prediction` | MSA-based (MMseqs2), diffusion decoder |
| Protenix | `protenix-prediction` | open-source AlphaFold3 reimplementation, MSA-based |
| OpenDDE | `opendde-prediction` | all-atom biomolecular complex prediction |
| ESMFold2 | `esmfold2-prediction` | all-atom, single-sequence (no MSA) |

All four ran on H100 in the Modal environment `proto-env`, stock tool configs, `seed=42`.

### Sequences

**Protein.** Mature chains only, taken from each UniProt entry's own `Chain` feature, so the
signal peptide (and for OSM the C-terminal propeptide) is stripped. This is what the models
should see: the signal peptide is cleaved before IL-6 ever encounters an aptamer, and leaving it
on would have the model fold 29 hydrophobic residues that do not exist in the mature cytokine.

| Target | UniProt | Full | Mature chain used | Length used |
|---|---|---|---|---|
| IL-6 | P05231 | 212 aa | 30–212 (signal 1–29 removed) | 183 aa |
| IL-11 | P20809 | 199 aa | 22–199 (signal 1–21 removed) | 178 aa |
| LIF | P15018 | 202 aa | 23–202 (signal 1–22 removed) | 180 aa |
| OSM | P13725 | 252 aa | 26–221 (signal 1–25 and propeptide 222–252 removed) | 196 aa |
| CNTF | P26441 | 200 aa | 1–200 (cytosolic, no signal peptide) | 200 aa |

**DNA.** Taken from the pipeline, never retyped: `parent45` from
`build_parents.parent("IL-6-7326.1")`, constructs from the `construct` column of
`switches.csv`. Because `switches.csv` is regenerated while the pipeline is being edited, it was
snapshotted at run start — `switches_snapshot.csv` and `inputs_snapshot.json` (with a SHA-256)
record exactly which version these predictions correspond to.

**Chain convention.** Chain `A` is always the protein; chain `B` is always the DNA.

## Files

| Path | Contents |
|---|---|
| `structures/<model>__<dna>__<target>.cif` | predicted complex, mmCIF |
| `results/<model>__<dna>__<target>.json` | one prediction: every metric, plus its contact epitope |
| `cofold_summary.json` | all predictions + cross-model agreement + off-target margins |
| `inputs_snapshot.json`, `switches_snapshot.csv` | exact inputs, with hash |
| `targets.json` | mature target sequences as fetched from UniProt |
| `run_cofold.py` | driver; resumable, one job at a time |
| `summarize.py` | builds `cofold_summary.json` from `results/` |

Reproduce with:

```bash
python run_cofold.py --phase A --models boltz2,protenix,opendde,esmfold2
python run_cofold.py --phase B --models boltz2
python summarize.py
```

## What the metrics mean

Recorded verbatim per prediction under `metrics`; not every model emits every field.

| Metric | Range | Better | Meaning |
|---|---|---|---|
| `iptm` | 0–1 | higher | **interface** predicted TM-score — the model's confidence in the *relative placement* of the two chains. This is the number that matters for "does the aptamer bind here", and the one to distrust most (see below). |
| `ptm` | 0–1 | higher | predicted TM-score for the whole complex, dominated by the protein's own fold |
| `chains_ptm` | 0–1 | higher | per-chain pTM; for chain `B` this is confidence in the DNA's *internal* fold |
| `pair_chains_iptm` | 0–1 | higher | ipTM per chain pair; the `A`/`B` off-diagonal entry is the protein–DNA interface |
| `complex_plddt` | 0–1 | higher | mean per-atom local confidence over the complex |
| `complex_iplddt` | 0–1 | higher | pLDDT restricted to interface atoms |
| `complex_pde` / `complex_ipde` | Å | lower | predicted distance error, all atoms / interface atoms |
| `avg_pae` / `pae` | 0–32 Å | lower | predicted aligned error; cross-chain blocks are the interface uncertainty |
| `confidence_score` | 0–1 | higher | each model's own blended headline score. Not comparable across models. |

**`epitope`** in each record is not a model output — it is measured from the returned
coordinates: every IL-6 residue with at least one heavy atom within **4.0 Å** of any DNA atom
(`Structure.interface_contact_residues`). This is the quantity compared across models, because
it is a geometric fact about the structure rather than a self-reported score.

**`cross_model_agreement`** reports, per DNA chain, the pairwise Jaccard overlap of those
epitope residue sets between models, plus the consensus set every model agrees on. High Jaccard
means independent models put the DNA in the same place on IL-6. That is the only confidence
statement here worth anything.

## Why these numbers do not rank anything

Predicted protein–DNA-aptamer structures are not currently accurate enough to rank binders, and
the published failure modes are specific and severe:

- **CASP16** included a 27-nt DNA aptamer target. **0 of 107** submitted models achieved under
  10 Å RMSD. Not "poor accuracy" — no group solved it.
- **AlphaFold3 memorises aptamers.** It averages **1.45 Å** on aptamer complexes published
  before its training cutoff and **6.40 Å** on ones published after. A good-looking prediction on
  a known aptamer is a lookup, not a generalisation, and every construct here is novel.
- **Interfaces are mostly wrong even when the fold is right.** Only **35.6%** of AF3
  protein–nucleic-acid predictions recover more than half the native contacts.

Consequences for reading this directory:

1. **ipTM is not an affinity and not a rank.** Do not order constructs by it. The pipeline's
   ranking comes from `switches.csv` (thermodynamics and measured Kd), and these predictions were
   run *after* that ranking, not into it.
2. **Cross-model agreement is the signal.** Four models that disagree about where the DNA sits
   tell you the binding mode is unresolved, regardless of how confident any one of them is. Four
   models that converge is weak positive evidence — they share training data and biases, so it is
   correlated evidence, not independent replication.
3. **The off-target margins are the most defensible output**, and only as a *relative*
   comparison. The same model, the same DNA, the same protocol, five homologous cytokines —
   systematic errors partly cancel, so the ordering is more trustworthy than any absolute value.
   A margin near zero means "this model cannot distinguish these targets", which is a real and
   useful negative result. It is not a measured specificity.
4. **The long poly-T linkers are not real structure.** Constructs 1–5 carry 27–54 nt of poly-T
   tether. Single-stranded poly-T is a flexible, largely unstructured spacer with no defined
   conformation, so whatever the models draw there is arbitrary — expect low pLDDT across the
   linker and ignore its geometry. Only the aptamer domain and the displaced strand are
   meaningful, and the tether's length can also drag down whole-complex scores like `ptm` and
   `complex_plddt` purely by adding disordered residues. Compare `parent45` (no tether) against
   the constructs with this in mind.
5. **Single-sequence vs MSA models are not equivalent.** ESMFold2 sees no MSA; Boltz-2, Protenix
   and OpenDDE do. Some of the between-model spread is that difference, not genuine uncertainty
   about the aptamer.
6. **No experimental validation is present in this directory.** Nothing here has been compared
   to a solved structure, a footprint, or a mutational scan. Use these as figures and as a
   hypothesis generator for which IL-6 surface to probe experimentally.

## Provenance caveats

- `switches.csv` was being regenerated by the pipeline during this run. The snapshot files pin
  what was actually used; re-running `summarize.py` against a newer `switches.csv` will not match.
- Modal initially refused to schedule H100/A100 on this workspace ("Please add a payment method
  to use H100 GPU functions"); only T4/L4/A10 were available. This was resolved mid-run and all
  reported predictions ran on stock `GPU_DEFAULT` (H100) with unmodified tool definitions.
