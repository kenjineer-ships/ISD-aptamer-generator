"""Aggregate results/*.json into cofold_summary.json + cross-model agreement stats.

    python summarize.py
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"

# Each tool's interface-confidence field, in preference order. Different models emit
# different names; we record all metrics verbatim and also normalise one comparable number.
IFACE_KEYS = ("iptm", "complex_iplddt", "interface_ptm", "ipae", "complex_ipde")


def iface_confidence(metrics: dict) -> tuple[str | None, float | None]:
    for k in IFACE_KEYS:
        v = metrics.get(k)
        if isinstance(v, (int, float)):
            return k, float(v)
    return None, None


def jaccard(a: set, b: set) -> float | None:
    if not a and not b:
        return None
    return len(a & b) / len(a | b)


TOL = 2


def jaccard_tolerant(a: set[int], b: set[int], tol: int = TOL) -> float | None:
    """Jaccard where a residue counts as shared if the other set has one within +/-tol.

    Residue-exact overlap is too strict for this comparison: two models can pick out the
    same surface patch and still score near zero because their contact lists are offset by
    a residue or two. This asks the weaker, more meaningful question - same site or not?
    """
    if not a and not b:
        return None
    a_hit = sum(any(x + d in b for d in range(-tol, tol + 1)) for x in a)
    b_hit = sum(any(x + d in a for d in range(-tol, tol + 1)) for x in b)
    denom = len(a) + len(b) - min(a_hit, b_hit)
    return (min(a_hit, b_hit) / denom) if denom else None


def segments(residues: set[int], gap: int = 4) -> list[str]:
    """Collapse a residue set into contiguous patches, e.g. ['11-18', '109-127']."""
    if not residues:
        return []
    xs = sorted(residues)
    out, start, prev = [], xs[0], xs[0]
    for x in xs[1:]:
        if x - prev > gap:
            out.append(f"{start}-{prev}" if start != prev else f"{start}")
            start = x
        prev = x
    out.append(f"{start}-{prev}" if start != prev else f"{start}")
    return out


def main() -> None:
    recs = []
    for p in sorted(RESULTS.glob("*.json")):
        recs.append(json.loads(p.read_text()))

    ok = [r for r in recs if r.get("success")]
    bad = [r for r in recs if not r.get("success")]

    for r in ok:
        k, v = iface_confidence(r.get("metrics", {}))
        r["interface_confidence"] = {"metric": k, "value": v}

    # ---- cross-model agreement (Phase A: all models, target IL6) ----
    phase_a = [r for r in ok if r.get("target") == "IL6" and r.get("phase") == "A"]
    epitopes: dict[tuple[str, str], set[int]] = {}
    for r in phase_a:
        ep = r.get("epitope", {})
        if "residues" in ep:
            epitopes[(r["dna_id"], r["model"])] = {int(x) for x in ep["residues"]}

    agreement = {"per_dna": {}, "note": "Jaccard overlap of the IL-6 residue set within 4.0 A of any DNA atom."}
    for dna_id in sorted({d for d, _ in epitopes}):
        models = sorted({m for d, m in epitopes if d == dna_id})
        pairs, pairs_tol = {}, {}
        for m1, m2 in combinations(models, 2):
            j = jaccard(epitopes[(dna_id, m1)], epitopes[(dna_id, m2)])
            pairs[f"{m1}|{m2}"] = None if j is None else round(j, 3)
            jt = jaccard_tolerant(epitopes[(dna_id, m1)], epitopes[(dna_id, m2)])
            pairs_tol[f"{m1}|{m2}"] = None if jt is None else round(jt, 3)
        vals = [v for v in pairs.values() if v is not None]
        vals_tol = [v for v in pairs_tol.values() if v is not None]
        # residues any model calls, vs residues every model calls
        sets = [epitopes[(dna_id, m)] for m in models]
        consensus = set.intersection(*sets) if sets else set()
        union = set.union(*sets) if sets else set()
        agreement["per_dna"][dna_id] = {
            "models": models,
            "n_contacts_per_model": {m: len(epitopes[(dna_id, m)]) for m in models},
            "contact_patches_per_model": {m: segments(epitopes[(dna_id, m)]) for m in models},
            "pairwise_jaccard": pairs,
            "mean_pairwise_jaccard": round(sum(vals) / len(vals), 3) if vals else None,
            "pairwise_jaccard_tolerant_pm2": pairs_tol,
            "mean_pairwise_jaccard_tolerant_pm2": round(sum(vals_tol) / len(vals_tol), 3) if vals_tol else None,
            "consensus_residues_all_models": sorted(consensus),
            "n_consensus": len(consensus),
            "n_union": len(union),
        }

    # Which model is the odd one out? Mean tolerant-Jaccard of each model against every
    # other model, averaged over DNA chains. A low score means that model puts the DNA
    # somewhere nobody else does.
    all_models = sorted({m for _, m in epitopes})
    outlier = {}
    for m in all_models:
        vals = []
        for dna_id in sorted({d for d, _ in epitopes}):
            for other in all_models:
                if other == m or (dna_id, other) not in epitopes or (dna_id, m) not in epitopes:
                    continue
                v = jaccard_tolerant(epitopes[(dna_id, m)], epitopes[(dna_id, other)])
                if v is not None:
                    vals.append(v)
        outlier[m] = round(sum(vals) / len(vals), 3) if vals else None
    agreement["mean_tolerant_jaccard_vs_all_others"] = outlier

    # Residues that recur as consensus across DNA chains, for the models that do agree.
    consensus_recurrence: dict[str, int] = {}
    for dna_id in sorted({d for d, _ in epitopes}):
        sets = [epitopes[(dna_id, m)] for m in all_models if (dna_id, m) in epitopes]
        if sets:
            for r in set.intersection(*sets):
                consensus_recurrence[str(r)] = consensus_recurrence.get(str(r), 0) + 1
    agreement["consensus_residue_recurrence_all_models"] = dict(
        sorted(consensus_recurrence.items(), key=lambda kv: -kv[1])
    )

    # ---- off-target margins (Phase B) ----
    phase_b = [r for r in ok if r.get("phase") == "B"]
    offtarget = {}
    by_model_dna: dict[tuple[str, str], list] = {}
    for r in phase_b:
        by_model_dna.setdefault((r["model"], r["dna_id"]), []).append(r)
    # on-target may live in Phase A; pull it in per (model, dna)
    for (model, dna_id), rows in by_model_dna.items():
        on = next((r for r in rows if r["target"] == "IL6"), None)
        if on is None:
            on = next(
                (r for r in ok if r["model"] == model and r["dna_id"] == dna_id and r["target"] == "IL6"),
                None,
            )
        if on is None:
            continue
        on_v = on["interface_confidence"]["value"]
        metric = on["interface_confidence"]["metric"]
        entry = {
            "model": model,
            "dna_id": dna_id,
            "metric": metric,
            "on_target_IL6": on_v,
            "off_targets": {},
        }
        for r in rows:
            if r["target"] == "IL6":
                continue
            v = r["interface_confidence"]["value"]
            entry["off_targets"][r["target"]] = {
                "uniprot": r["target_uniprot"],
                "value": v,
                "margin_vs_IL6": None if (v is None or on_v is None) else round(on_v - v, 4),
                "n_contact_residues": r.get("epitope", {}).get("n_contact_residues"),
            }
        entry["min_margin"] = min(
            (d["margin_vs_IL6"] for d in entry["off_targets"].values() if d["margin_vs_IL6"] is not None),
            default=None,
        )
        offtarget[f"{model}__{dna_id}"] = entry

    summary = {
        "interpretation": (
            "These co-folds are FIGURES AND A WEAK SPECIFICITY SIGNAL, NOT RANKING INPUT. The "
            "aptamer pipeline deliberately does not rank on predicted structure. CASP16's 27-nt "
            "DNA aptamer target had 0 of 107 submitted models under 10 A RMSD; AF3 averages 1.45 A "
            "on aptamers published before its training cutoff vs 6.40 A after it; only 35.6% of AF3 "
            "protein-nucleic-acid predictions recover more than half the native contacts. Treat "
            "cross-model agreement as the signal, never a single model's own confidence score."
        ),
        "counts": {
            "total_predictions_attempted": len(recs),
            "succeeded": len(ok),
            "failed": len(bad),
            "phase_A": len([r for r in recs if r.get("phase") == "A"]),
            "phase_B": len([r for r in recs if r.get("phase") == "B"]),
        },
        "chain_convention": {"A": "protein target (mature chain)", "B": "DNA aptamer / switch construct"},
        "cross_model_agreement": agreement,
        "off_target_specificity": offtarget,
        "failures": [
            {k: r.get(k) for k in ("key", "model", "dna_id", "target", "error")} for r in bad
        ],
        "predictions": recs,
    }
    (HERE / "cofold_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"{len(ok)} ok / {len(bad)} failed -> cofold_summary.json")
    for dna_id, a in agreement["per_dna"].items():
        print(
            f"  {dna_id}: meanJaccard={a['mean_pairwise_jaccard']} "
            f"tolerant+-2={a['mean_pairwise_jaccard_tolerant_pm2']} "
            f"consensus={a['n_consensus']}/{a['n_union']}"
        )
        for m, segs in a["contact_patches_per_model"].items():
            print(f"      {m:9s} {','.join(segs)}")
    print(f"  outlier check (mean tol-Jaccard vs others): {agreement['mean_tolerant_jaccard_vs_all_others']}")
    for k, e in offtarget.items():
        print(f"  {k} [{e['metric']}] on-target={e['on_target_IL6']} min_margin={e['min_margin']}")
        for t, d in e["off_targets"].items():
            print(f"      {t}: {d['value']} (margin {d['margin_vs_IL6']})")


if __name__ == "__main__":
    main()
