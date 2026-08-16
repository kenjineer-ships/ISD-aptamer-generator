"""Co-fold aptamer/switch DNA chains against IL-6 and gp130-family off-targets.

Runs one (model, dna, target) job at a time, writes each result to its own JSON under
results/ and each structure under structures/, and skips jobs already done. Safe to
re-run: it resumes.

    python run_cofold.py --phase A --models boltz2,protenix,opendde,esmfold2
    python run_cofold.py --phase B --models boltz2

These predictions are FIGURES AND A WEAK SPECIFICITY SIGNAL, NOT RANKING INPUT.
See README.md in this directory.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
APTAMER = HERE.parent
RESULTS = HERE / "results"
STRUCTS = HERE / "structures"

# UniProt accession -> (label, mature-chain start, mature-chain end) as annotated by
# UniProt's own Chain feature. 1-indexed, inclusive.
TARGETS = {
    "P05231": ("IL6", 30, 212),  # signal peptide 1-29
    "P20809": ("IL11", 22, 199),  # signal peptide 1-21
    "P15018": ("LIF", 23, 202),  # signal peptide 1-22
    "P13725": ("OSM", 26, 221),  # signal 1-25, C-terminal propeptide 222-252 removed
    "P26441": ("CNTF", 1, 200),  # cytosolic, no signal peptide
}

MODELS = {
    "boltz2": ("proto_tools.tools.structure_prediction.boltz2.boltz2", "run_boltz2", "Boltz2Input", "Boltz2Config"),
    "protenix": (
        "proto_tools.tools.structure_prediction.protenix.protenix",
        "run_protenix",
        "ProtenixInput",
        "ProtenixConfig",
    ),
    "opendde": (
        "proto_tools.tools.structure_prediction.opendde.opendde",
        "run_opendde",
        "OpenDDEInput",
        "OpenDDEConfig",
    ),
    "esmfold2": (
        "proto_tools.tools.structure_prediction.esmfold2.esmfold2",
        "run_esmfold2",
        "ESMFold2Input",
        "ESMFold2Config",
    ),
}

PROTEIN_CHAIN = "A"
DNA_CHAIN = "B"
SEED = 42
CONTACT_CUTOFF = 4.0


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def fetch_targets() -> dict[str, dict]:
    """Mature target sequences, cached to targets.json."""
    cache = HERE / "targets.json"
    if cache.exists():
        return json.loads(cache.read_text())
    out = {}
    for acc, (label, start, end) in TARGETS.items():
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        with urllib.request.urlopen(url, timeout=60) as r:
            lines = r.read().decode().splitlines()
        full = "".join(x.strip() for x in lines if not x.startswith(">"))
        mature = full[start - 1 : end]
        out[label] = {
            "label": label,
            "uniprot": acc,
            "full_len": len(full),
            "mature_range": [start, end],
            "mature_seq": mature,
            "mature_len": len(mature),
            "note": f"UniProt Chain feature {start}-{end} of {acc}; signal peptide/propeptide stripped",
        }
        log(f"{label} ({acc}): full {len(full)} aa -> mature {len(mature)} aa ({start}-{end})")
    cache.write_text(json.dumps(out, indent=2))
    return out


def load_dna_chains() -> list[dict]:
    """Bare 45-nt parent + top 5 ranked constructs, from the snapshot taken at run start."""
    snap = json.loads((HERE / "inputs_snapshot.json").read_text())
    chains = [
        {
            "dna_id": "parent45",
            "sequence": snap["parent"]["sequence"],
            "kind": "bare parent aptamer",
            "rank": 0,
            "meta": {"kd_M": snap["parent"]["kd_M"], "name": snap["parent"]["name"]},
        }
    ]
    for row in snap["top5"]:
        chains.append(
            {
                "dna_id": f"construct{row['rank']}",
                "sequence": row["construct"],
                "kind": "switch construct",
                "rank": row["rank"],
                "meta": {k: row[k] for k in row if k != "construct"},
            }
        )
    return chains


def job_key(model: str, dna_id: str, target: str) -> str:
    return f"{model}__{dna_id}__{target}"


def run_one(model: str, dna: dict, target: dict) -> dict:
    """Run a single co-fold. Returns a result record."""
    import importlib

    mod_path, run_name, in_name, cfg_name = MODELS[model]
    mod = importlib.import_module(mod_path)
    run_fn = getattr(mod, run_name)
    InputCls = getattr(mod, in_name)
    CfgCls = getattr(mod, cfg_name)

    complex_spec = {
        "chains": [
            {"id": PROTEIN_CHAIN, "sequence": target["mature_seq"], "entity_type": "protein"},
            {"id": DNA_CHAIN, "sequence": dna["sequence"], "entity_type": "dna"},
        ]
    }

    # device="modal" dispatches to the app deployed in the Modal workspace. Without it the
    # call runs locally, which on Windows fails outright ("Unsupported operating system").
    cfg_kwargs = {"seed": SEED, "verbose": 1, "device": "modal"}
    # Not every model exposes every knob; only pass what its config declares.
    for k, v in (("timeout", 3600),):
        if k in CfgCls.model_fields:
            cfg_kwargs[k] = v
    cfg = CfgCls(**{k: v for k, v in cfg_kwargs.items() if k in CfgCls.model_fields})

    key = job_key(model, dna["dna_id"], target["label"])
    log(f"RUN {key}  protein={target['mature_len']}aa dna={len(dna['sequence'])}nt")
    t0 = time.time()
    out = run_fn(InputCls(complexes=[complex_spec]), cfg)
    wall = time.time() - t0

    st = out.structures[0]
    cif_path = STRUCTS / f"{key}.cif"
    st.write_cif(cif_path)

    metrics = {}
    if st.metrics is not None:
        raw = st.metrics.model_dump() if hasattr(st.metrics, "model_dump") else dict(st.metrics)
        # pae is O(n^2) and not requested; drop it from the summary if it showed up.
        metrics = {k: v for k, v in raw.items() if k != "pae"}

    rec = {
        "key": key,
        "phase": "A" if target["label"] == "IL6" else "B",
        "model": model,
        "dna_id": dna["dna_id"],
        "dna_chain_id": DNA_CHAIN,
        "dna_sequence": dna["sequence"],
        "dna_len": len(dna["sequence"]),
        "dna_kind": dna["kind"],
        "dna_rank": dna["rank"],
        "dna_meta": dna["meta"],
        "target": target["label"],
        "target_uniprot": target["uniprot"],
        "target_chain_id": PROTEIN_CHAIN,
        "target_seq_len": target["mature_len"],
        "target_mature_range": target["mature_range"],
        "seed": SEED,
        "wall_seconds": round(wall, 1),
        "success": True,
        "structure_file": f"structures/{key}.cif",
        "metrics": metrics,
    }

    # Where does the DNA touch the protein? This is the cross-model agreement signal.
    try:
        contacts = st.interface_contact_residues(PROTEIN_CHAIN, DNA_CHAIN, cutoff=CONTACT_CUTOFF)
        rec["epitope"] = {
            "cutoff_angstrom": CONTACT_CUTOFF,
            "n_contact_residues": len(contacts),
            "residues": {str(k): v for k, v in sorted(contacts.items())},
        }
    except Exception as e:  # analysis must never lose a finished structure
        rec["epitope"] = {"error": f"{type(e).__name__}: {e}"}

    try:
        plddt = st.per_residue_plddt
        if plddt:
            vals = list(plddt.values()) if isinstance(plddt, dict) else list(plddt)
            nums = [float(v) for v in vals if isinstance(v, (int, float))]
            if nums:
                rec["mean_per_residue_plddt"] = round(sum(nums) / len(nums), 4)
    except Exception:
        pass

    try:
        rec["chain_ids_in_output"] = list(st.get_chain_ids())
    except Exception:
        pass

    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=["A", "B"], required=True)
    ap.add_argument("--models", required=True, help="comma-separated model slugs")
    ap.add_argument("--best-construct", default=None, help="phase B: dna_id to use (default: construct1)")
    args = ap.parse_args()

    RESULTS.mkdir(parents=True, exist_ok=True)
    STRUCTS.mkdir(parents=True, exist_ok=True)

    targets = fetch_targets()
    dnas = load_dna_chains()
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    for m in models:
        if m not in MODELS:
            log(f"unknown model {m}")
            return 2

    jobs: list[tuple[str, dict, dict]] = []
    if args.phase == "A":
        for m in models:
            for d in dnas:
                jobs.append((m, d, targets["IL6"]))
    else:
        best = args.best_construct or "construct1"
        d = next((x for x in dnas if x["dna_id"] == best), None)
        if d is None:
            log(f"no such dna_id {best}")
            return 2
        for m in models:
            for label in ("IL6", "IL11", "LIF", "OSM", "CNTF"):
                jobs.append((m, d, targets[label]))

    log(f"phase {args.phase}: {len(jobs)} job(s) over models {models}")
    done = failed = skipped = 0
    for model, dna, target in jobs:
        key = job_key(model, dna["dna_id"], target["label"])
        path = RESULTS / f"{key}.json"
        if path.exists():
            log(f"SKIP {key} (already done)")
            skipped += 1
            continue
        try:
            rec = run_one(model, dna, target)
            path.write_text(json.dumps(rec, indent=2))
            ip = rec["metrics"].get("iptm")
            ne = rec.get("epitope", {}).get("n_contact_residues")
            log(f"OK   {key}  iptm={ip} contacts={ne} {rec['wall_seconds']}s")
            done += 1
        except Exception as e:
            tb = traceback.format_exc()
            path.write_text(
                json.dumps(
                    {
                        "key": key,
                        "phase": args.phase,
                        "model": model,
                        "dna_id": dna["dna_id"],
                        "dna_sequence": dna["sequence"],
                        "target": target["label"],
                        "target_uniprot": target["uniprot"],
                        "success": False,
                        "error": f"{type(e).__name__}: {e}",
                        "traceback": tb[-4000:],
                    },
                    indent=2,
                )
            )
            log(f"FAIL {key}: {type(e).__name__}: {e}")
            failed += 1

    log(f"phase {args.phase} finished: {done} ok, {failed} failed, {skipped} skipped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
