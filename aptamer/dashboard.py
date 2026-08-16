"""Render the pipeline's outputs as one self-contained dashboard.html at the repo root.

Reads whatever the pipeline has written so far (parents.json, switches.csv, mismatches.csv,
negative_controls.csv, *.pdb / *.cif) and inlines it into a single HTML file as JSON. Nothing
is fetched at runtime: the dashboard is opened over file://, where Chrome blocks local
XHR/fetch as cross-origin.

    python aptamer/dashboard.py        -> ./dashboard.html

Note which parent the constructs sit on. switch_library.PARENT is "IL-6-7326.1", a 45-nt
TRUNCATION (module B removed) of the 74-nt "IL-6-7326" row in parents.json -- that name
matches no row of the parent table, so construct_parent() resolves it through
build_parents.parent() and the page says so in four places (page banner, the parent table's
role column with the removed tail struck through, the switch-library heading, views 9-10).
A reader must never be left to assume the 74-mer.

Missing inputs degrade to a "not yet generated" card; they never raise. No science is
recomputed here -- dg_switch, KD_app, structures and MFEs are read as produced. The only
derived values are display conveniences (residence time 1/koff, fractional occupancy from a
published K_D), and those come from build_parents so the formulas live in one place.

Charts are inline SVG built by ~200 lines of vanilla JS from the inlined JSON. Two panels
prefer a CDN script for something SVG-by-hand cannot do: 3Dmol.js for the 3D structures
(view 8) and fornac for the secondary-structure layout (view 7). View 7 degrades to the
built-in base-pair arc diagram when fornac is unreachable, so only view 8 actually needs
internet; everything else works with no network at all.

Views 11-12 read aptamer/cofold/ (written by run_cofold.py + summarize.py, and never written
here -- that directory is another step's output). Both are negative results and are drawn as
such: the four models do not agree on where the aptamer binds, and Boltz-2 ipTM ranks CNTF
above IL-6 for the one chain with a measured K_D. Only five of the 32 .cif complexes are
inlined into view 8; see COFOLD_STRUCTURES for which five and why.

fornac ships its stylesheet inside the bundle and style-loader injects it into <head>
unscoped, where its bare `svg {}` / `text {}` rules would resize every other chart. A
scoped copy is inlined from vendor/fornac.css by scoped_fornac_css() and the injected
copy is deleted at runtime -- see dropFornacGlobalCss() in the JS.
"""

import csv
import json
import pathlib
import statistics
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "dashboard.html"

# --- pipeline constants: import them, never re-type them ----------------------------------
sys.path.insert(0, str(HERE))

IL6_MW = 23718  # fallback only; overwritten by the import below
WOUND = {
    "post-surgical (median)": 986.6,
    "inflammatory (median)": 4964,
    "infected (median)": 5883,
    "max observed": 135500,
}
SKIN_MODEL_PG_ML = 5000.0  # 5 ng/mL, README "Known constraint: the affinity gap"
NOTES = []


def _occupancy_fallback(kd, pg_ml):
    conc = pg_ml * 1e-9 / IL6_MW
    return conc, conc / (kd + conc)


occupancy = _occupancy_fallback
try:
    from build_parents import IL6_MW as _MW, WOUND as _WOUND, occupancy as _occ

    IL6_MW, WOUND, occupancy = _MW, _WOUND, _occ
except Exception as exc:  # ViennaRNA missing, etc. -- display values still work
    NOTES.append(
        f"build_parents not importable ({exc.__class__.__name__}: {exc}); using local constants")

SWITCH_PARENT = "IL-6-7326"
RANDOMISED = sorted({12, 13, 18, 19, 20, 25, 26, 27, 45, 46, 47, 53, 54, 55, 60, 61})
# The design target as switch_library states it. Older revisions exposed it as a literal
# DG_WINDOW; newer ones state the affinity budget it is derived from. Take whichever exists
# and never re-derive it here -- the shaded band falls back to the passing data range.
DG_WINDOW = None
BUDGET = {}
try:
    import switch_library as _sl

    SWITCH_PARENT = _sl.PARENT
    RANDOMISED = sorted(_sl.randomised_positions())
    if hasattr(_sl, "DG_WINDOW"):
        DG_WINDOW = tuple(_sl.DG_WINDOW)
    for k in ("KD_APP_MAX_NM", "CLOSED_MIN", "MIN_ENGAGEMENT"):
        if hasattr(_sl, k):
            BUDGET[k] = getattr(_sl, k)
except Exception as exc:
    NOTES.append(
        f"switch_library not importable ({exc.__class__.__name__}: {exc}); "
        "parent name and randomised block positions taken from the cached values"
    )


# --- readers -----------------------------------------------------------------------------
def read_parents():
    p = HERE / "parents.json"
    if not p.exists():
        return None
    try:
        rows = json.loads(p.read_text(encoding="utf-8"))
    except Exception as exc:
        NOTES.append(f"parents.json unreadable: {exc}")
        return None
    for r in rows:
        koff = r.get("koff_s") or 0
        r["residence_s"] = (1.0 / koff) if koff else None
        r["kd_nM"] = (r.get("KD_M") or 0) * 1e9
    return rows


NUMERIC = {
    "ds_len": int, "linker_len": int, "rand_covered": int,
    "gc": float, "dg_switch": float, "closed_frac": float,
    "kd_app_nM": float, "engagement": float,
}
# Columns newer revisions of switch_library added. Cast them when present but never require
# them: a switches.csv written before they existed must still fill the table, so a missing
# key here costs one feature (the view-9 tether comparison) rather than every switch row.
OPTIONAL = {
    "selectivity": float, "tether_nt": int,
    "dg_designed": float, "dg_homodimer": float, "dg_hairpin": float,
}


def read_switches():
    p = HERE / "switches.csv"
    if not p.exists():
        return None
    try:
        with p.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
    except Exception as exc:
        NOTES.append(f"switches.csv unreadable: {exc}")
        return None
    out = []
    for i, r in enumerate(rows):
        try:
            for k, cast in NUMERIC.items():
                r[k] = cast(r[k])
            for k, cast in OPTIONAL.items():
                if r.get(k) not in (None, ""):
                    try:
                        r[k] = cast(r[k])
                    except Exception:
                        r.pop(k, None)
            start, end = r["window"].split("-")
            r["w0"], r["w1"] = int(start), int(end)
            r["rank"] = i + 1  # file order IS the pipeline ranking
            out.append(r)
        except Exception:
            continue  # a malformed row is dropped, not fatal
    if len(out) != len(rows):
        NOTES.append(f"switches.csv: {len(rows) - len(out)} malformed row(s) skipped")
    return out


def _read_rows(name, numeric, script):
    """CSV -> list of dicts with `numeric` cast. Missing file returns None; a malformed row is
    dropped and reported. Nothing here raises: the section it feeds degrades instead."""
    p = HERE / name
    if not p.exists():
        return None
    try:
        with p.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
    except Exception as exc:
        NOTES.append(f"{name} unreadable: {exc}")
        return None
    out = []
    for i, r in enumerate(rows):
        try:
            for k, cast in numeric.items():
                r[k] = cast(r[k])
            r["rank"] = i + 1  # file order IS the pipeline's own ordering
            out.append(r)
        except Exception:
            continue
    if len(out) != len(rows):
        NOTES.append(f"{name}: {len(rows) - len(out)} malformed row(s) skipped -- "
                     f"regenerate with python aptamer/{script}")
    return out


MM_NUMERIC = {
    "linker_len": int, "mismatch_pos": int, "rand_covered": int, "tether_nt": int,
    "kd_app_wt_nM": float, "d_kd_nM": float, "dg_switch": float, "closed_frac": float,
    "kd_app_nM": float, "engagement": float, "selectivity": float,
}
NC_NUMERIC = {
    "linker_len": int, "tether_nt": int, "dg_switch": float, "closed_frac": float,
    "kd_app_nM": float, "engagement": float, "selectivity": float,
    "passes": lambda v: str(v).strip().lower() in ("true", "1", "yes"),
}
ARMS = ("designed", "scrambled", "reversed", "foreign")
ARM_NOTE = {
    "designed": "the displacement strand the pipeline designed",
    "scrambled": "same bases, shuffled — complementarity destroyed",
    "reversed": "same bases, reversed not complemented — wrong pairing register",
    "foreign": "reverse complement of a window of a SHUFFLED aptamer — a real duplex-former "
               "that is simply not complementary to this aptamer (the strict control)",
}


def read_mismatches():
    return _read_rows("mismatches.csv", MM_NUMERIC, "mismatch_tune.py")


def read_controls():
    return _read_rows("negative_controls.csv", NC_NUMERIC, "negative_controls.py")


def mismatch_context(mismatches, switches):
    """Price each variant's affinity in tether length instead of mismatches.

    A mismatch is the only knob that moves K_D,app without moving either length, so the honest
    comparison is: what tether would the unmutated library need to reach the same K_D,app while
    covering at least as many randomised nucleotides? Annotates each row with `tether_equiv`
    and returns the headline numbers. None when switches.csv predates tether_nt.
    """
    if not mismatches:
        return None
    pool = [r for r in (switches or [])
            if isinstance(r.get("tether_nt"), int) and "rand_covered" in r]
    if not pool:
        NOTES.append("switches.csv has no tether_nt column, so view 9 shows the before/after "
                     "without the equivalent-tether comparison")
        return None
    for r in mismatches:
        ts = [q["tether_nt"] for q in pool
              if q["kd_app_nM"] <= r["kd_app_nM"] and q["rand_covered"] >= r["rand_covered"]]
        r["tether_equiv"] = min(ts) if ts else None
    tethers = [r["tether_nt"] for r in mismatches]
    covered = max(r["rand_covered"] for r in mismatches)
    best = min(mismatches, key=lambda r: r["kd_app_nM"])
    best_cov = min((r for r in mismatches if r["rand_covered"] == covered),
                   key=lambda r: r["kd_app_nM"])
    lib = [q for q in pool if q["rand_covered"] >= covered]
    lib_best = min(lib, key=lambda q: q["kd_app_nM"]) if lib else None
    # The construct the mismatch is actually competing with: the SHORTEST-tether unmutated
    # construct that matches best_cov's affinity while covering at least as much. lib_best is a
    # different (and weaker) comparison -- it is the tightest anywhere in the library, and it
    # gets there by spending tether.
    alt = [q for q in lib if q["kd_app_nM"] <= best_cov["kd_app_nM"]]
    match_alt = min(alt, key=lambda q: (q["tether_nt"], q["kd_app_nM"])) if alt else None
    keep = ("ds", "window", "kd_app_nM", "tether_nt", "rand_covered")
    return {
        "n": len(mismatches),
        "nBase": len({(r["ds_wt"], r["linker_len"]) for r in mismatches}),
        "improved": sum(1 for r in mismatches if r["d_kd_nM"] < 0),
        "tetherLo": min(tethers), "tetherHi": max(tethers),
        "best": {k: best[k] for k in
                 ("name", "kd_app_wt_nM", "kd_app_nM", "d_kd_nM", "tether_nt",
                  "rand_covered", "tether_equiv")},
        "covered": covered,
        "bestCovered": {k: best_cov[k] for k in
                        ("name", "kd_app_wt_nM", "kd_app_nM", "d_kd_nM", "tether_nt",
                         "rand_covered", "tether_equiv")},
        "libBest": None if not lib_best else {k: lib_best[k] for k in keep},
        "matchAlt": None if not match_alt else {k: match_alt[k] for k in keep},
    }


def control_summary(controls):
    """Median score per arm plus pass rate. Medians, not means: the control arms are floored at
    engagement 0 and a mean would be dragged by the handful of accidental pairings."""
    if not controls:
        return []
    arms = [a for a in ARMS if any(r["arm"] == a for r in controls)]
    arms += sorted({r["arm"] for r in controls} - set(arms))
    med = lambda rows, k: statistics.median(r[k] for r in rows)  # noqa: E731
    designed = [r for r in controls if r["arm"] == "designed"]
    ref = med(designed, "engagement") if designed else None
    out = []
    for a in arms:
        rows = [r for r in controls if r["arm"] == a]
        out.append({
            "arm": a, "n": len(rows), "note": ARM_NOTE.get(a, ""),
            "dg_switch": round(med(rows, "dg_switch"), 2),
            "closed_frac": round(med(rows, "closed_frac"), 3),
            "kd_app_nM": round(med(rows, "kd_app_nM"), 1),
            "engagement": round(med(rows, "engagement"), 2),
            "engMax": round(max(r["engagement"] for r in rows), 2),
            "dgLo": round(min(r["dg_switch"] for r in rows), 2),
            "dgHi": round(max(r["dg_switch"] for r in rows), 2),
            "passRate": sum(1 for r in rows if r["passes"]) / len(rows),
            "reachRef": None if ref is None else
            sum(1 for r in rows if r["engagement"] >= ref) / len(rows),
        })
    return out


def construct_parent(parents, switches):
    """The sequence the constructs are actually built on -- resolved, never assumed.

    switch_library.PARENT is "IL-6-7326.1", which is not a row in parents.json: it is the
    45-nt truncation of "IL-6-7326". Preference order: build_parents.parent(), the single
    place that owns the truncation; else peel a construct apart (construct = aptamer + linker
    + DS, so the aptamer is a straight prefix slice); else the exact-name row if one exists.
    """
    seq = None
    try:
        import build_parents as _bp

        seq = _bp.parent(SWITCH_PARENT)[0]
    except Exception:
        pass
    if not seq:
        for r in switches or []:
            n = len(r.get("construct", "")) - r["linker_len"] - r["ds_len"]
            if n > 0:
                seq = r["construct"][:n]
                break
    if not seq:
        seq = next((p["sequence"] for p in (parents or []) if p["name"] == SWITCH_PARENT), None)
    if not seq:
        NOTES.append(f"could not resolve the sequence of the construct parent {SWITCH_PARENT}; "
                     "view 6 and the parent labelling degrade")
        return None
    exact = next((p for p in (parents or []) if p["sequence"] == seq), None)
    src = next((p for p in (parents or [])
                if p["sequence"].startswith(seq) and p["length"] > len(seq)), None)
    return {
        "name": SWITCH_PARENT,
        "sequence": seq,
        "length": len(seq),
        "truncated": bool(src),
        "source": src["name"] if src else (exact["name"] if exact else None),
        "sourceLength": src["length"] if src else None,
        "removed": src["sequence"][len(seq):] if src else "",
        "kdNM": (src or exact or {}).get("kd_nM"),
    }


# --- co-folded complexes: aptamer/cofold/, READ ONLY ---------------------------------------
# That directory belongs to the co-folding step (run_cofold.py + summarize.py). Nothing here
# writes to it and nothing here re-derives its science: ipTM, the 4.0 A contact-residue sets
# and the Jaccard matrices are read exactly as summarize.py produced them.
COFOLD_DIR = HERE / "cofold"
# Lane order for view 11. boltz2 / esmfold2 / protenix broadly converge and opendde does not,
# so the trio sits together and the outlier sits last: the ordering is part of the figure.
COFOLD_MODELS = ("boltz2", "esmfold2", "protenix", "opendde")
TRIO = ("boltz2", "esmfold2", "protenix")
# The complexes inlined into view 8. Each .cif is 200-250 KB and 32 of them sit on disk, so
# inlining the lot would add ~8 MB to the page. These five are the ones that carry an argument:
# one DNA chain (parent45 -- the only chain with a measured K_D) placed by all four models, so
# the disagreement is visible in 3D, plus the CNTF run that outscores the true target. The
# other 27 stay in aptamer/cofold/structures/ and the panel says so.
COFOLD_STRUCTURES = (
    "boltz2__parent45__IL6", "esmfold2__parent45__IL6",
    "protenix__parent45__IL6", "opendde__parent45__IL6",
    "boltz2__parent45__CNTF",
)
TARGET_LABEL = {"IL6": "IL-6", "IL11": "IL-11", "LIF": "LIF", "OSM": "OSM", "CNTF": "CNTF"}


def _f(v):
    """float(v) or None -- dna_meta values arrive as strings straight out of the CSV."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def read_cofold():
    """cofold_summary.json -> everything views 11 and 12 draw.

    Two display aggregations are computed here and nothing else: how many of the four models
    call each residue (a tally of the contact sets already in the file, the same shape as the
    view-6 coverage count), and the boltz2/esmfold2/protenix intersection that the cofold
    README names as the recurring patch. Every number that is a measurement -- ipTM, pTM,
    pLDDT, the residue sets, the Jaccard matrices, the off-target margins -- is read.
    """
    p = COFOLD_DIR / "cofold_summary.json"
    if not p.exists():
        return None
    try:
        s = json.loads(p.read_text(encoding="utf-8"))
    except Exception as exc:
        NOTES.append(f"cofold/cofold_summary.json unreadable ({exc.__class__.__name__}: {exc}); "
                     "views 11-12 show 'not yet generated'")
        return None
    preds = [r for r in (s.get("predictions") or []) if r.get("success")]
    if not preds:
        NOTES.append("cofold_summary.json has no successful predictions; views 11-12 degrade")
        return None
    by_key = {r.get("key"): r for r in preds}
    agree = s.get("cross_model_agreement") or {}
    per = agree.get("per_dna") or {}

    def ordered(models):
        return ([m for m in COFOLD_MODELS if m in models]
                + sorted(set(models) - set(COFOLD_MODELS)))

    chains, trio_rec, target_len, offset = [], {}, None, None
    # parent45 leads: it is the only chain with a measured K_D, so it is the one that can fail
    # a control. The constructs follow in pipeline rank order.
    for name in sorted(per, key=lambda n: (0 if n.startswith("parent") else 1, n)):
        a = per[name] or {}
        entries = []
        for m in ordered(a.get("models") or []):
            rec = by_key.get(f"{m}__{name}__IL6")
            if not rec:
                continue
            ep = (rec.get("epitope") or {}).get("residues") or {}
            met = rec.get("metrics") or {}
            cptm = met.get("chains_ptm") or []
            target_len = target_len or rec.get("target_seq_len")
            rng = rec.get("target_mature_range") or []
            if offset is None and len(rng) == 2:
                offset = int(rng[0]) - 1  # mature residue 1 == UniProt residue rng[0]
            entries.append({
                "model": m,
                "residues": sorted(int(k) for k in ep),
                "aa": {str(int(k)): v for k, v in ep.items()},
                "n": (a.get("n_contacts_per_model") or {}).get(m, len(ep)),
                "patches": (a.get("contact_patches_per_model") or {}).get(m) or [],
                "iptm": met.get("iptm"), "ptm": met.get("ptm"),
                "dnaPtm": cptm[-1] if cptm else None,
                "plddt": met.get("complex_plddt"), "conf": met.get("confidence_score"),
            })
        if not entries:
            continue
        counts = {}
        for e in entries:
            for r in e["residues"]:
                counts[r] = counts.get(r, 0) + 1
        trio_sets = [set(e["residues"]) for e in entries if e["model"] in TRIO]
        trio = sorted(set.intersection(*trio_sets)) if len(trio_sets) == len(TRIO) else []
        for r in trio:
            trio_rec[r] = trio_rec.get(r, 0) + 1
        rec0 = by_key.get(f"{entries[0]['model']}__{name}__IL6") or {}
        meta = rec0.get("dna_meta") or {}
        chains.append({
            "dna": name,
            "len": rec0.get("dna_len"),
            "kind": rec0.get("dna_kind"),
            "rank": rec0.get("dna_rank"),
            "kdAppNM": _f(meta.get("kd_app_nM")),
            "tether": _f(meta.get("tether_nt")),
            "models": entries,
            "counts": {str(k): v for k, v in sorted(counts.items())},
            "trio": trio,
            "jaccard": a.get("pairwise_jaccard") or {},
            "jaccardTol": a.get("pairwise_jaccard_tolerant_pm2") or {},
            "meanJ": a.get("mean_pairwise_jaccard"),
            "meanJTol": a.get("mean_pairwise_jaccard_tolerant_pm2"),
            "consensus": [int(r) for r in (a.get("consensus_residues_all_models") or [])],
            "nConsensus": a.get("n_consensus"),
            "nUnion": a.get("n_union"),
        })
    if not chains:
        return None

    # --- Phase B: the off-target panels, one per DNA chain --------------------------------
    off = []
    for key, o in (s.get("off_target_specificity") or {}).items():
        dna, model = o.get("dna_id"), o.get("model")
        on = o.get("on_target_IL6")
        onrec = by_key.get(f"{model}__{dna}__IL6") or {}
        rows = [{
            "target": "IL6", "label": TARGET_LABEL["IL6"],
            "uniprot": onrec.get("target_uniprot") or "P05231",
            "iptm": on, "margin": 0.0, "onTarget": True,
            "nContacts": (onrec.get("epitope") or {}).get("n_contact_residues"),
        }]
        for t, v in sorted((o.get("off_targets") or {}).items()):
            rows.append({
                "target": t, "label": TARGET_LABEL.get(t, t), "uniprot": v.get("uniprot"),
                "iptm": v.get("value"), "margin": v.get("margin_vs_IL6"), "onTarget": False,
                "nContacts": v.get("n_contact_residues"),
            })
        rows.sort(key=lambda r: -(r["iptm"] or 0))
        off.append({
            "key": key, "dna": dna, "model": model, "metric": o.get("metric") or "iptm",
            "onTarget": on, "rows": rows, "minMargin": o.get("min_margin"),
            "dnaLen": onrec.get("dna_len"), "dnaKind": onrec.get("dna_kind"),
            # a name here is the metric failing its control, not a specificity result
            "beats": [r["label"] for r in rows
                      if not r["onTarget"] and (r["iptm"] or 0) > (on or 0)],
        })
    off.sort(key=lambda d: (0 if str(d["dna"]).startswith("parent") else 1, d["dna"]))

    sdir = COFOLD_DIR / "structures"
    return {
        "counts": s.get("counts") or {},
        "interpretation": s.get("interpretation"),
        "contactNote": agree.get("note"),
        "chains": chains,
        "models": ordered({e["model"] for c in chains for e in c["models"]}),
        "trioModels": [m for m in TRIO],
        "vsOthers": agree.get("mean_tolerant_jaccard_vs_all_others") or {},
        "allModelRecurrence": agree.get("consensus_residue_recurrence_all_models") or {},
        "trioRecurrence": sorted(trio_rec.items(), key=lambda kv: (-kv[1], kv[0])),
        "targetLen": target_len,
        "uniprotOffset": 29 if offset is None else offset,
        "offTarget": off,
        "nOnDisk": len(list(sdir.glob("*.cif"))) if sdir.is_dir() else 0,
    }


def read_structures(cofold=None):
    """1ALU, anything the structure step dropped in aptamer/, plus the five inlined co-folds."""
    out = []
    for p in sorted(HERE.glob("*.pdb")) + sorted(HERE.glob("*.cif")):
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except Exception as exc:
            NOTES.append(f"{p.name} unreadable: {exc}")
            continue
        is_ref = p.stem.upper() == "1ALU"
        out.append({
            "name": p.name,
            "format": "pdb" if p.suffix.lower() == ".pdb" else "cif",
            "kind": "receptor" if is_ref else "complex",
            "label": "human IL-6 crystal structure, 1.9 A (PDB 1ALU)" if is_ref else p.stem,
            "text": text,
            "atoms": sum(1 for line in text.splitlines() if line.startswith(("ATOM", "HETATM"))),
            "chainA": None, "chainB": None,
        })
    return out + read_cofold_structures(cofold)


def read_cofold_structures(cofold=None):
    """The five .cif complexes from COFOLD_STRUCTURES, labelled from their own run records.

    Chain A is the protein and chain B the DNA in every file the co-folding step writes, which
    is what lets view 8 colour them apart. A missing file is skipped, not fatal.
    """
    d = COFOLD_DIR / "structures"
    if not d.is_dir():
        return []
    meta, dna_len = {}, {}
    for c in (cofold or {}).get("chains", []):
        dna_len[c["dna"]] = c.get("len")
        for e in c["models"]:
            meta[f"{e['model']}__{c['dna']}__IL6"] = (c, e["iptm"])
    # the off-target runs use the same DNA chains, so their length is known even though the
    # per-chain agreement record only covers the on-target runs
    for panel in (cofold or {}).get("offTarget", []):
        for r in panel["rows"]:
            meta.setdefault(f"{panel['model']}__{panel['dna']}__{r['target']}",
                            ({"len": dna_len.get(panel["dna"])}, r["iptm"]))
    out, missing = [], []
    for stem in COFOLD_STRUCTURES:
        p = d / (stem + ".cif")
        if not p.exists():
            missing.append(stem)
            continue
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except Exception as exc:
            NOTES.append(f"cofold/structures/{p.name} unreadable: {exc}")
            continue
        try:
            model, dna, target = stem.split("__")
        except ValueError:
            model, dna, target = stem, "?", "?"
        chain, iptm = meta.get(stem, (None, None))
        nt = f", {chain['len']} nt" if chain and chain.get("len") else ""
        tag = " — ipTM %.3f" % iptm if isinstance(iptm, (int, float)) else ""
        out.append({
            "name": p.name,
            "format": "cif",
            "kind": "complex",
            "label": (f"{model} · {dna}{nt} vs {TARGET_LABEL.get(target, target)}"
                      f"{tag}{' (off-target control)' if target != 'IL6' else ''}"),
            "text": text,
            "atoms": sum(1 for line in text.splitlines() if line.startswith(("ATOM", "HETATM"))),
            "chainA": TARGET_LABEL.get(target, target) + " (protein)",
            "chainB": f"{dna} DNA{nt}",
        })
    if missing:
        NOTES.append("cofold/structures: missing " + ", ".join(m + ".cif" for m in missing)
                     + " -- view 8 lists the ones that are there")
    return out


def occupancy_series(parents):
    """Fraction bound vs IL-6 (pg/mL), log-spaced. Display derivation only."""
    if not parents:
        return []
    series = []
    for p in parents:
        kd = p.get("KD_M")
        if not kd:
            continue
        pts = []
        for i in range(121):
            pg = 10.0 ** (1 + 6 * i / 120)  # 1e1 .. 1e7 pg/mL
            _, frac = occupancy(kd, pg)
            pts.append([pg, frac])
        half = kd * IL6_MW / 1e-9  # pg/mL giving 50% occupancy == K_D in mass units
        series.append({"name": p["name"], "kd_nM": p["kd_nM"], "pts": pts, "half_pg": half})
    return series


def occupancy_table(parents):
    rows = []
    landmarks = list(WOUND.items()) + [("skin model (5 ng/mL)", SKIN_MODEL_PG_ML)]
    for label, pg in landmarks:
        row = {"label": label, "pg_ml": pg, "pM": pg * 1e-9 / IL6_MW * 1e12, "frac": {}}
        for p in parents or []:
            if p.get("KD_M"):
                row["frac"][p["name"]] = occupancy(p["KD_M"], pg)[1]
        rows.append(row)
    return rows


def design_window(switches):
    """(lo, hi, caption) for the shaded ΔG_switch band. Stated target if the pipeline
    exposes one, otherwise the range that actually passed -- never re-derived here."""
    if DG_WINDOW:
        return DG_WINDOW[0], DG_WINDOW[1], (
            f"design window {DG_WINDOW[0]} to {DG_WINDOW[1]} kcal/mol")
    if not switches:
        return -2.2, -0.7, "design window (not stated by the pipeline)"
    dgs = [r["dg_switch"] for r in switches]
    bits = []
    if "KD_APP_MAX_NM" in BUDGET:
        bits.append(f"K_D,app ≤ {BUDGET['KD_APP_MAX_NM']:.0f} nM")
    if "CLOSED_MIN" in BUDGET:
        bits.append(f"closed ≥ {BUDGET['CLOSED_MIN']:.0%}")
    cap = "ΔG range that passed" + (" (" + ", ".join(bits) + ")" if bits else "")
    return min(dgs), max(dgs), cap


def cofold_staleness(cofold, switches, mismatches, cparent):
    """How far the pipeline has moved since the co-folds were run.

    The five constructs in aptamer/cofold/ were snapshotted before the ensemble-ΔG and
    K_closed corrections and before mismatch tuning existed, so their lengths and K_D,app are
    not the shortlist's any more. This puts the two side by side from live files rather than
    from a number typed into the caveat -- if the pipeline moves again, the caveat moves with
    it. None when there is nothing to compare against.
    """
    if not cofold:
        return None
    used = [c for c in cofold["chains"] if c.get("kdAppNM")]
    if not used:
        return None
    top = min(used, key=lambda c: c.get("rank") or 99)
    # a mismatch variant carries no `construct` column, so its length is the parent plus the
    # tether the row records -- the same arithmetic switch_library does when it builds one
    base = (cparent or {}).get("length")
    cand = [{"kd": r["kd_app_nM"], "len": len(r.get("construct") or "") or None,
             "what": f"switch construct DS {r['ds']} at window {r['window']}"}
            for r in switches or []]
    for r in mismatches or []:
        teth = r.get("tether_nt")
        cand.append({"kd": r["kd_app_nM"],
                     "len": (base + teth) if (base and isinstance(teth, int)) else None,
                     "what": f"single-mismatch variant {r.get('name', '')}".strip()})
    now = min(cand, key=lambda c: c["kd"]) if cand else None
    # The length half of the comparison: the shortest thing the pipeline now has that already
    # beats every co-folded chain on K_D,app. Affinity alone understates the move -- the
    # constructs got much shorter as well as tighter.
    ceiling = min(c["kdAppNM"] for c in used)
    beat = [c for c in cand if c["kd"] < ceiling and c["len"]]
    short = min(beat, key=lambda c: (c["len"], c["kd"])) if beat else None
    return {
        "short": short,
        "then": {"dna": top["dna"], "len": top.get("len"), "kd": top.get("kdAppNM"),
                 "rank": top.get("rank")},
        "thenKdLo": min(c["kdAppNM"] for c in used),
        "thenKdHi": max(c["kdAppNM"] for c in used),
        "thenLenLo": min(c["len"] for c in used if c.get("len")),
        "thenLenHi": max(c["len"] for c in used if c.get("len")),
        "nUsed": len(used),
        "now": now,
    }


def build_data():
    parents = read_parents()
    switches = read_switches()
    mismatches = read_mismatches()
    controls = read_controls()
    cofold = read_cofold()
    structures = read_structures(cofold)
    lo, hi, cap = design_window(switches)
    cparent = construct_parent(parents, switches)
    return {
        "generated": None,
        "parents": parents,
        "switches": switches,
        "mismatches": mismatches,
        "mismatchHeadline": mismatch_context(mismatches, switches),
        "controls": controls,
        "controlSummary": control_summary(controls),  # each row carries its own arm note
        "structures": structures,
        "cofold": cofold,
        # what a construct in the co-fold set cost then vs what the pipeline reaches now --
        # views 11-12 predate the ensemble-dG / K_closed corrections and the mismatch tuning
        "cofoldStale": cofold_staleness(cofold, switches, mismatches, cparent),
        "occupancy": occupancy_series(parents),
        "occupancyTable": occupancy_table(parents),
        "wound": WOUND,
        "skinModel": SKIN_MODEL_PG_ML,
        "il6MW": IL6_MW,
        "dgWindow": [lo, hi],
        "dgWindowCaption": cap,
        "budget": BUDGET,
        "randomised": RANDOMISED,
        "switchParent": SWITCH_PARENT,
        # The truncated parent, resolved -- an exact-name lookup in parents.json returns nothing
        # for "IL-6-7326.1" and used to leave view 6 empty.
        "constructParent": cparent,
        "parentSeq": cparent["sequence"] if cparent else None,
        "notes": NOTES,
    }


# --- the page ----------------------------------------------------------------------------
CSS = r"""
:root { color-scheme: light; }
.viz-root {
  --surface-1:#fcfcfb; --page:#f9f9f7;
  --text-primary:#0b0b0b; --text-secondary:#52514e; --muted:#898781;
  --grid:#e1e0d9; --axis:#c3c2b7; --border:rgba(11,11,11,0.10);
  --series-1:#2a78d6; --series-2:#eb6834; --series-3:#1baf7a; --series-4:#eda100;
  --warning:#fab219; --critical:#d03b3b; --good:#0ca30c;
  --wash:rgba(42,120,214,0.10);
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) .viz-root {
    color-scheme: dark;
    --surface-1:#1a1a19; --page:#0d0d0d;
    --text-primary:#ffffff; --text-secondary:#c3c2b7; --muted:#898781;
    --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,0.10);
    --series-1:#3987e5; --series-2:#d95926; --series-3:#199e70; --series-4:#c98500;
    --wash:rgba(57,135,229,0.14);
  }
}
:root[data-theme="dark"] .viz-root {
  color-scheme: dark;
  --surface-1:#1a1a19; --page:#0d0d0d;
  --text-primary:#ffffff; --text-secondary:#c3c2b7; --muted:#898781;
  --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,0.10);
  --series-1:#3987e5; --series-2:#d95926; --series-3:#199e70; --series-4:#c98500;
  --wash:rgba(57,135,229,0.14);
}

* { box-sizing:border-box; }
html, body { margin:0; }
body { background:var(--page); }
.viz-root {
  font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
  background:var(--page); color:var(--text-primary);
  padding:24px clamp(16px,4vw,48px) 96px; min-height:100vh;
}
h1 { font-size:22px; font-weight:600; margin:0; letter-spacing:-0.01em; }
h2 { font-size:15px; font-weight:600; margin:0 0 2px; }
p.sub, .sub { color:var(--text-secondary); font-size:13px; margin:0; }
.muted { color:var(--muted); }
a { color:var(--series-1); }

header.top { display:flex; align-items:flex-start; gap:16px; flex-wrap:wrap;
  padding-bottom:16px; border-bottom:1px solid var(--border); margin-bottom:20px; }
header.top .grow { flex:1 1 320px; }
button, select, input[type=number], input[type=text] {
  font:inherit; font-size:13px; color:var(--text-primary); background:var(--surface-1);
  border:1px solid var(--border); border-radius:8px; padding:6px 10px;
}
button { cursor:pointer; }
button:hover, select:hover { border-color:var(--axis); }

.card { background:var(--surface-1); border:1px solid var(--border); border-radius:12px;
  padding:16px 18px 18px; margin:0 0 20px; }
.card > header { margin-bottom:12px; }
.grid2 { display:grid; grid-template-columns:repeat(auto-fit,minmax(420px,1fr)); gap:20px; }
.stats { display:grid; grid-template-columns:repeat(auto-fit,minmax(180px,1fr)); gap:12px;
  margin-bottom:20px; }
.stat { background:var(--surface-1); border:1px solid var(--border); border-radius:12px;
  padding:14px 16px; }
.stat .label { font-size:12px; color:var(--text-secondary); }
.stat .value { font-size:26px; font-weight:600; line-height:1.2; margin-top:2px; }
.stat .value.hero { font-size:48px; }
.stat .unit { font-size:13px; color:var(--text-secondary); font-weight:400; margin-left:4px; }

table { border-collapse:separate; border-spacing:0; width:100%; font-size:13px;
  font-variant-numeric:tabular-nums; }
th, td { text-align:right; padding:7px 10px; white-space:nowrap;
  border-bottom:1px solid var(--grid); }
th:first-child, td:first-child, th.l, td.l { text-align:left; }
thead th { position:sticky; top:0; background:var(--surface-1); z-index:2;
  color:var(--text-secondary); font-weight:600; border-bottom:1px solid var(--axis);
  font-size:12px; }
thead th.sortable { cursor:pointer; user-select:none; }
thead th.sortable:hover { color:var(--text-primary); }
thead th .arrow { color:var(--muted); font-size:10px; }
tbody tr:hover { background:var(--wash); }
.scroll { max-height:460px; overflow:auto; border:1px solid var(--border);
  border-radius:10px; }
code, .mono { font-family: ui-monospace, "Cascadia Mono", Consolas, monospace; font-size:12px; }
.seq { font-family: ui-monospace, "Cascadia Mono", Consolas, monospace; font-size:11.5px;
  letter-spacing:0.055em; word-break:break-all; }

.badge { display:inline-flex; align-items:center; gap:5px; font-size:11.5px; font-weight:600;
  border-radius:999px; padding:2px 9px 2px 7px; border:1px solid var(--border);
  color:var(--text-primary); }
.badge.warn { background:rgba(250,178,25,0.18); }
.badge.ok { background:rgba(12,163,12,0.14); color:var(--text-secondary); font-weight:500; }
.badge .ic { font-weight:700; }
.issues { margin:10px 0 0; padding:10px 12px; border-radius:10px;
  background:rgba(250,178,25,0.14); border:1px solid rgba(250,178,25,0.45); font-size:12.5px; }
.issues b { font-weight:600; }
.issues ul { margin:6px 0 0 18px; padding:0; }
/* caveats that must be read before the chart, not after it: same amber as .issues, but sits
   above the figure it qualifies. Views 9-10 use it; the text is static HTML so it survives
   even if the chart JS fails. */
.caveat { margin:0 0 12px; padding:10px 12px; border-radius:10px;
  background:rgba(250,178,25,0.14); border:1px solid rgba(250,178,25,0.45); font-size:12.5px;
  line-height:1.5; }
.caveat + .caveat { margin-top:-4px; }
.caveat b { font-weight:600; }
.callout { margin:0 0 14px; padding:11px 13px; border-radius:10px; background:var(--wash);
  border:1px solid var(--border); font-size:12.5px; line-height:1.55; }
.callout b { font-weight:600; }
/* a struck-through nucleotide run: the part of a full-length parent the construct drops */
.cut { color:var(--muted); text-decoration:line-through; text-decoration-thickness:1px; }
.keep { background:var(--wash); border-radius:3px; }
/* a colour key inline in a table cell -- the number itself stays in an ink token */
.dot { width:8px; height:8px; border-radius:50%; display:inline-block; margin-right:6px;
  vertical-align:middle; }

.filters { display:flex; flex-wrap:wrap; gap:14px 18px; align-items:flex-end;
  background:var(--surface-1); border:1px solid var(--border); border-radius:12px;
  padding:12px 16px; margin-bottom:16px; }
.filters .f { display:flex; flex-direction:column; gap:4px; }
.filters label { font-size:11.5px; color:var(--text-secondary); }
.filters .row { display:flex; gap:6px; align-items:center; }
.filters input[type=number] { width:88px; }
.filters .count { font-size:12.5px; color:var(--text-secondary); margin-left:auto; }

.legend { display:flex; flex-wrap:wrap; gap:6px 16px; font-size:12px;
  color:var(--text-secondary); margin:8px 0 2px; }
.legend .k { display:inline-flex; align-items:center; gap:6px; }
.legend .sw { width:12px; height:12px; border-radius:3px; display:inline-block; }
.legend .ln { width:16px; height:2px; border-radius:2px; display:inline-block; }
.rampbar { display:flex; align-items:center; gap:8px; font-size:12px;
  color:var(--text-secondary); }
.rampbar .bar { height:10px; width:160px; border-radius:5px; border:1px solid var(--border); }

svg { display:block; max-width:100%; overflow:visible; }
svg text { fill:var(--text-secondary); }
.tip { position:fixed; z-index:50; pointer-events:none; opacity:0; transition:opacity .08s;
  background:var(--surface-1); color:var(--text-primary); border:1px solid var(--border);
  box-shadow:0 6px 24px rgba(0,0,0,.16); border-radius:10px; padding:8px 10px;
  font-size:12.5px; max-width:320px; }
.tip .v { font-weight:600; font-variant-numeric:tabular-nums; }
.tip .n { color:var(--text-secondary); }
.tip .r { display:flex; align-items:center; gap:7px; }
.tip .k { width:14px; height:2px; border-radius:2px; flex:0 0 auto; }
.empty { border:1px dashed var(--axis); border-radius:10px; padding:22px;
  color:var(--text-secondary); font-size:13px; text-align:center; }
details.tv { margin-top:10px; }
details.tv summary { font-size:12.5px; color:var(--text-secondary); cursor:pointer; }
.viewer { position:relative; height:420px; border-radius:10px; overflow:hidden;
  border:1px solid var(--border); background:var(--surface-1); }
footer.foot { color:var(--muted); font-size:12px; border-top:1px solid var(--border);
  padding-top:14px; margin-top:8px; }
@media print { .viewer { display:none; } }
"""

# --- fornac's stylesheet, rewritten so it can only ever touch view 7 ----------------------
FORNA_SCOPE = "#ss-forna"  # the div fornaView() creates inside section 7's #ss-chart


def scoped_fornac_css(scope=FORNA_SCOPE):
    """Read vendor/fornac.css and prefix every selector with `scope`.

    The vendored sheet contains bare `svg { width:100%; min-height:100% }` and
    `text { pointer-events:none }` rules. Every other chart on this page is hand-built
    inline SVG, so an unscoped copy would visually destroy views 2, 4, 5 and 6. The file is
    flat (no at-rules, no nesting), so a split on braces is a sufficient parser -- anything
    unexpected is dropped and reported as a generator note rather than passed through.
    """
    p = HERE / "vendor" / "fornac.css"
    try:
        raw = p.read_text(encoding="utf-8")
    except Exception as exc:
        NOTES.append(f"vendor/fornac.css unreadable ({exc.__class__.__name__}: {exc}); "
                     f"fornac is not loaded and view 7 falls back to the arc diagram")
        return ""
    blocks, dropped = [], []
    for chunk in raw.split("}"):
        if "{" not in chunk:
            continue  # trailing whitespace after the last rule
        sel, body = chunk.split("{", 1)
        sel, body = " ".join(sel.split()), " ".join(body.split())
        if not sel or sel.startswith("@") or "{" in body:
            dropped.append(sel or "(empty selector)")
            continue
        parts = [s.strip() for s in sel.split(",") if s.strip()]
        if not parts:
            continue
        blocks.append(", ".join(f"{scope} {s}" for s in parts) + " { " + body + " }")
    if dropped:
        NOTES.append("vendor/fornac.css: could not scope " + ", ".join(dropped) +
                     " -- rule(s) dropped")
    # our own overrides, after the vendored rules so they win at equal specificity
    blocks += [
        f"{scope} {{ position:relative; overflow:hidden; border-radius:10px;",
        "  border:1px solid var(--border); background:var(--surface-1); }",
        # fornac.css sets `color` on the nucleotide labels, which does nothing to SVG text,
        # so the dashboard's own `svg text` rule would win. fornac tags each label with
        # label_type, which lets the letters and the position numbers be styled apart.
        # Every region colour in ssRegions() is dark enough for white letters (inkOn()).
        f"{scope} svg {{ overflow:hidden; }}",
        f"{scope} text {{ fill:var(--text-secondary); }}",
        f"{scope} text.node-label[label_type=nucleotide] {{ fill:#ffffff; }}",
        f"{scope} text.node-label[label_type=label] {{ fill:var(--text-secondary); }}",
        f"{scope} circle.node.label {{ fill:var(--surface-1); }}",
        f"{scope} circle.node {{ stroke:var(--border); }}",
        # fornac ships raw red / blue / green links; pull them into the page palette.
        # `line.link` is the backbone; base pairs carry both classes, and come later here,
        # so they win. Neither hue collides with the region colours on the nucleotides.
        f"{scope} line.link {{ stroke:var(--muted); }}",
        f"{scope} line.basepair {{ stroke:var(--text-primary); stroke-width:2; }}",
        f"{scope} line.pseudoknot {{ stroke:var(--series-4); }}",
        f"{scope} line.intermolecule {{ stroke:var(--series-3); }}",
    ]
    return "\n".join(blocks)


JS = r"""
'use strict';
const DATA = JSON.parse(document.getElementById('dash-data').textContent);
const NS = 'http://www.w3.org/2000/svg';

/* ---------- palette: both modes are selected, not flipped ---------- */
const PAL = {
  light: { s1:'#2a78d6', s2:'#eb6834', s3:'#1baf7a', s4:'#eda100',
           surface:'#fcfcfb', grid:'#e1e0d9', axis:'#c3c2b7', ink:'#0b0b0b',
           ink2:'#52514e', muted:'#898781', warn:'#fab219', crit:'#d03b3b',
           /* sequential blue 100->700 */
           seq:['#cde2fb','#b7d3f6','#9ec5f4','#86b6ef','#6da7ec','#5598e7','#3987e5',
                '#2a78d6','#256abf','#1c5cab','#184f95','#104281','#0d366b'],
           /* ordinal blue: nothing lighter than step 250 on the light surface */
           ord:['#86b6ef','#6da7ec','#3987e5','#2a78d6','#256abf','#1c5cab','#104281'] },
  dark:  { s1:'#3987e5', s2:'#d95926', s3:'#199e70', s4:'#c98500',
           surface:'#1a1a19', grid:'#2c2c2a', axis:'#383835', ink:'#ffffff',
           ink2:'#c3c2b7', muted:'#898781', warn:'#fab219', crit:'#d03b3b',
           seq:['#0d366b','#104281','#184f95','#1c5cab','#256abf','#2a78d6','#3987e5',
                '#5598e7','#6da7ec','#86b6ef','#9ec5f4','#b7d3f6','#cde2fb'],
           /* ordinal blue on dark: nothing darker than step 600 */
           ord:['#b7d3f6','#86b6ef','#5598e7','#3987e5','#2a78d6','#256abf','#184f95'] }
};
function mode() {
  const t = document.documentElement.dataset.theme;
  if (t === 'dark' || t === 'light') return t;
  return matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
}
let P = PAL[mode()];
const SERIES = m => [m.s1, m.s2, m.s3, m.s4];
function seqColor(t) {          /* t in [0,1] -> sequential blue */
  const r = P.seq, i = Math.max(0, Math.min(r.length - 1, Math.round(t * (r.length - 1))));
  return r[i];
}
function ordColor(i, n) {       /* ordered classes -> ordinal blue */
  const r = P.ord, k = n <= 1 ? 0 : Math.round(i / (n - 1) * (r.length - 1));
  return r[k];
}
function inkOn(hex) {           /* label inside a fill: pick by luminance */
  const h = hex.replace('#',''), n = parseInt(h, 16);
  const l = (0.2126*((n>>16)&255) + 0.7152*((n>>8)&255) + 0.0722*(n&255)) / 255;
  return l > 0.55 ? '#0b0b0b' : '#ffffff';
}

/* ---------- tiny SVG + scale helpers ---------- */
function el(tag, attrs, parent) {
  const n = document.createElementNS(NS, tag);
  for (const k in (attrs || {})) if (attrs[k] !== null && attrs[k] !== undefined)
    n.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(n);
  return n;
}
function txt(parent, x, y, s, attrs) {
  const t = el('text', Object.assign({x:x, y:y, 'font-size':11}, attrs || {}), parent);
  t.appendChild(document.createTextNode(s));
  return t;
}
const lin = (d0,d1,r0,r1) => v => r0 + (v - d0) / (d1 - d0) * (r1 - r0);
const log = (d0,d1,r0,r1) => {
  const a = Math.log10(d0), b = Math.log10(d1);
  return v => r0 + (Math.log10(Math.max(v, 1e-300)) - a) / (b - a) * (r1 - r0);
};
const fmt = (v, d) => v.toLocaleString(undefined, {minimumFractionDigits:d===undefined?0:d,
                                                  maximumFractionDigits:d===undefined?1:d});
const sup = n => ({'-':'⁻','0':'⁰','1':'¹','2':'²','3':'³','4':'⁴',
                   '5':'⁵','6':'⁶','7':'⁷','8':'⁸','9':'⁹'})[n] || n;
const pow10 = e => '10' + String(e).split('').map(sup).join('');

/* ---------- one tooltip for the whole page ---------- */
const TIP = document.createElement('div');
TIP.className = 'tip'; document.body.appendChild(TIP);
function showTip(ev, rows) {
  TIP.textContent = '';
  rows.forEach(r => {
    const d = document.createElement('div');
    d.className = 'r';
    if (r.color) { const k = document.createElement('span');
      k.className = 'k'; k.style.background = r.color; d.appendChild(k); }
    if (r.value !== undefined) { const v = document.createElement('span');
      v.className = 'v'; v.textContent = r.value; d.appendChild(v); }
    if (r.name) { const n = document.createElement('span');
      n.className = 'n'; n.textContent = r.name; d.appendChild(n); }
    TIP.appendChild(d);
  });
  TIP.style.opacity = 1;
  const pad = 14, w = TIP.offsetWidth, h = TIP.offsetHeight;
  let x = ev.clientX + pad, y = ev.clientY + pad;
  if (x + w > innerWidth - 8) x = ev.clientX - w - pad;
  if (y + h > innerHeight - 8) y = ev.clientY - h - pad;
  TIP.style.left = x + 'px'; TIP.style.top = y + 'px';
}
const hideTip = () => { TIP.style.opacity = 0; };
addEventListener('scroll', hideTip, {passive:true});

function grid(g, x1, x2, y, first) {
  el('line', {x1:x1, x2:x2, y1:y, y2:y, stroke:first ? P.axis : P.grid,
              'stroke-width':1, 'shape-rendering':'crispEdges'}, g);
}
function card(id) { return document.getElementById(id); }
function empty(node, msg) {
  node.textContent = '';
  const d = document.createElement('div');
  d.className = 'empty'; d.textContent = msg; node.appendChild(d);
}
function guard(id, fn) {
  const n = card(id); if (!n) return;
  try { fn(n); } catch (e) { empty(n, 'could not render this view: ' + e.message); }
}

/* ---------- a sortable table in ~25 lines; no datatable library -------------------------
   cols: {k, h, d (decimals), text (left-align + monospace), cell (custom renderer), tip}
   state: {k, d, draw} -- `draw` re-renders the caller so the arrow and order stay in sync. */
function sortableTable(node, cols, rows, state) {
  const num = v => (v === null || v === undefined) ? Infinity : v;  /* blanks sort last */
  const data = rows.slice().sort((a, b) => {
    const x = a[state.k], y = b[state.k];
    if (typeof x === 'string' || typeof y === 'string')
      return state.d * String(x).localeCompare(String(y));
    return state.d * (num(x) - num(y));
  });
  node.textContent = '';
  const wrap = document.createElement('div'); wrap.className = 'scroll';
  const t = document.createElement('table');
  const thead = document.createElement('thead'), tr = document.createElement('tr');
  cols.forEach(c => {
    const th = document.createElement('th');
    th.className = 'sortable' + (c.text ? ' l' : '');
    th.textContent = c.h;
    const a = document.createElement('span'); a.className = 'arrow';
    a.textContent = state.k === c.k ? (state.d > 0 ? ' ▲' : ' ▼') : '';
    th.appendChild(a);
    th.title = c.tip || ('sort by ' + c.h);
    th.addEventListener('click', () => {
      if (state.k === c.k) state.d = -state.d; else { state.k = c.k; state.d = 1; }
      state.draw();
    });
    tr.appendChild(th);
  });
  thead.appendChild(tr); t.appendChild(thead);
  const tb = document.createElement('tbody');
  data.forEach(r => {
    const row = document.createElement('tr');
    cols.forEach(c => {
      const td = document.createElement('td');
      if (c.cell) { c.cell(td, r); }
      else if (c.text) { td.className = 'l mono'; td.textContent = r[c.k]; }
      else td.textContent = (r[c.k] === null || r[c.k] === undefined) ? '—' : fmt(r[c.k], c.d);
      row.appendChild(td);
    });
    if (r.__title) row.title = r.__title;
    tb.appendChild(row);
  });
  t.appendChild(tb); wrap.appendChild(t); node.appendChild(wrap);
  if (!data.length) empty(node, 'no rows');
}

/* =========================================================================
   VIEW 1 - parent aptamer table
   ========================================================================= */
function renderParents() {
  guard('parents-body', node => {
    const rows = DATA.parents;
    if (!rows || !rows.length) return empty(node, 'parents.json not yet generated - run python aptamer/build_parents.py');
    node.textContent = '';
    const wrap = document.createElement('div'); wrap.className = 'scroll';
    const t = document.createElement('table');
    t.innerHTML = '';
    const head = ['Aptamer', 'Variant', 'Length (nt)', 'K_D (nM)', 'k_off (s⁻¹)',
                  'Residence time 1/k_off (s)', 'MFE 22 °C (kcal/mol)',
                  'MFE 37 °C (kcal/mol)', 'Reconstruction',
                  'Role in this pipeline'];
    const thead = document.createElement('thead'), htr = document.createElement('tr');
    head.forEach((h, i) => { const th = document.createElement('th');
      th.textContent = h; if (i < 2 || i >= 8) th.className = 'l'; htr.appendChild(th); });
    thead.appendChild(htr); t.appendChild(thead);
    const tb = document.createElement('tbody');
    const cp = DATA.constructParent;
    rows.forEach(r => {
      const isSource = !!(cp && cp.truncated && r.name === cp.source);
      const isParent = !!(cp && r.name === cp.name);
      const tr = document.createElement('tr');
      const cells = [
        [r.name, 'l'], [r.variant, 'l'], [String(r.length), ''],
        [fmt(r.kd_nM, 1), ''], [(r.koff_s || 0).toExponential(2), ''],
        [r.residence_s ? fmt(r.residence_s, 0) : '—', ''],
        [(r.mfe_22C).toFixed(2), ''], [(r.mfe_37C).toFixed(2), ''],
        ['@badge', 'l'], ['@role', 'l']
      ];
      cells.forEach(([v, cls]) => {
        const td = document.createElement('td'); if (cls) td.className = cls;
        if (v === '@badge') {
          const b = document.createElement('span');
          const bad = (r.issues || []).length;
          b.className = 'badge ' + (bad ? 'warn' : 'ok');
          const ic = document.createElement('span'); ic.className = 'ic';
          ic.textContent = bad ? '⚠' : '✓'; b.appendChild(ic);
          b.appendChild(document.createTextNode(bad ? bad + ' uncertain nt' : 'clean'));
          td.appendChild(b);
        } else if (v === '@role') {
          /* the point of this column: the constructs are built on a 45-nt truncation, so no
             row of this table is itself the scaffold. Say which row it came from. */
          if (isSource) {
            const b = document.createElement('span'); b.className = 'badge warn';
            const ic = document.createElement('span'); ic.className = 'ic';
            ic.textContent = '⚠'; b.appendChild(ic);
            b.appendChild(document.createTextNode('full-length source, not the scaffold'));
            td.appendChild(b);
            const s = document.createElement('div');
            s.className = 'sub'; s.style.marginTop = '3px';
            s.textContent = 'constructs use its first ' + cp.length + ' nt as ' + cp.name;
            td.appendChild(s);
          } else if (isParent) {
            td.textContent = 'the construct parent (' + r.length + ' nt)';
          } else {
            td.className = 'l muted';
            td.textContent = cp ? 'not used in the switch library' : '—';
          }
        } else td.textContent = v;
        tr.appendChild(td);
      });
      tb.appendChild(tr);
      const sr = document.createElement('tr');
      const td = document.createElement('td');
      td.colSpan = head.length; td.className = 'l';
      const s = document.createElement('div'); s.className = 'seq';
      if (isSource) {
        /* kept vs removed, drawn on the sequence itself: the strike-through is the module the
           truncation drops, so the 45-mer the constructs use is visible in place. */
        const keep = document.createElement('span');
        keep.className = 'keep'; keep.textContent = r.sequence.slice(0, cp.length);
        const drop = document.createElement('span');
        drop.className = 'cut'; drop.textContent = r.sequence.slice(cp.length);
        s.appendChild(keep); s.appendChild(drop);
      } else s.textContent = r.sequence;
      td.appendChild(s);
      const d = document.createElement('div');
      d.className = 'seq muted'; d.style.marginTop = '2px';
      d.textContent = r.structure_37C + '   (37 °C)';
      td.appendChild(d);
      if (isSource) {
        const c = document.createElement('div');
        c.className = 'sub'; c.style.marginTop = '4px';
        c.textContent = 'kept 1–' + cp.length + ' (highlighted) → ' + cp.name +
          ' · removed ' + (cp.length + 1) + '–' + r.length + ' (struck through, module B)';
        td.appendChild(c);
      }
      sr.appendChild(td); tb.appendChild(sr);
    });
    t.appendChild(tb); wrap.appendChild(t); node.appendChild(wrap);

    const bad = rows.filter(r => (r.issues || []).length);
    if (bad.length) {
      const box = document.createElement('div'); box.className = 'issues';
      const b = document.createElement('b');
      b.textContent = '⚠ Do not order these blind — reconstruction warnings';
      box.appendChild(b);
      const ul = document.createElement('ul');
      bad.forEach(r => (r.issues || []).forEach(i => {
        const li = document.createElement('li');
        li.textContent = r.name + ': ' + i;
        ul.appendChild(li);
      }));
      box.appendChild(ul);
      node.appendChild(box);
    }
  });
}

/* =========================================================================
   VIEW 2 - occupancy curve
   ========================================================================= */
function renderOccupancy() {
  guard('occ-chart', node => {
    const S = DATA.occupancy;
    if (!S || !S.length) return empty(node, 'parents.json not yet generated — no occupancy curve');
    node.textContent = '';
    const W = Math.max(520, Math.min(node.clientWidth || 760, 980)), H = 380;
    const m = {t:14, r:120, b:46, l:52};
    const svg = el('svg', {width:'100%', viewBox:'0 0 ' + W + ' ' + H,
                           role:'img', 'aria-label':'Fraction of aptamer bound vs IL-6 concentration'}, node);
    const x = log(10, 1e7, m.l, W - m.r), y = lin(0, 1, H - m.b, m.t);

    /* wound-fluid bands: context washes, drawn under the data */
    const bands = Object.entries(DATA.wound);
    bands.forEach(([label, pg], i) => {
      const px = x(pg);
      el('line', {x1:px, x2:px, y1:m.t, y2:H - m.b, stroke:P.axis, 'stroke-width':1,
                  'shape-rendering':'crispEdges'}, svg);
      const lab = txt(svg, px + 4, m.t + 12 + (i % 2) * 13, label.replace(' (median)', ''),
                      {'font-size':10.5, fill:P.muted});
      lab.setAttribute('text-anchor', 'start');
    });
    el('rect', {x:x(bands[0][1]), y:m.t, width:Math.max(1, x(bands[2][1]) - x(bands[0][1])),
                height:H - m.b - m.t, fill:P.s1, opacity:0.07}, svg);
    /* skin-model reading */
    const sx = x(DATA.skinModel);
    el('line', {x1:sx, x2:sx, y1:m.t, y2:H - m.b, stroke:P.warn, 'stroke-width':2}, svg);
    txt(svg, sx + 5, H - m.b - 6, 'skin model 5 ng/mL', {'font-size':10.5, fill:P.ink2});

    /* axes */
    for (let f = 0; f <= 1.0001; f += 0.25) {
      grid(svg, m.l, W - m.r, y(f), f === 0);
      txt(svg, m.l - 8, y(f) + 4, Math.round(f * 100) + '%',
          {'text-anchor':'end', 'font-size':10.5, fill:P.muted});
    }
    for (let e = 1; e <= 7; e++) {
      txt(svg, x(Math.pow(10, e)), H - m.b + 16, pow10(e),
          {'text-anchor':'middle', 'font-size':10.5, fill:P.muted});
    }
    txt(svg, (m.l + W - m.r) / 2, H - 8, 'IL-6 (pg/mL, log scale)',
        {'text-anchor':'middle', 'font-size':11.5, fill:P.ink2});
    txt(svg, m.l - 8, m.t - 2, 'fraction bound', {'text-anchor':'end', 'font-size':11.5, fill:P.ink2});

    /* lines + a direct label at each curve's 50% crossing (its K_D) */
    const cols = SERIES(P);
    S.forEach((s, i) => {
      const d = s.pts.map((p, k) => (k ? 'L' : 'M') + x(p[0]).toFixed(1) + ' ' + y(p[1]).toFixed(1)).join(' ');
      el('path', {d:d, fill:'none', stroke:cols[i % 4], 'stroke-width':2,
                  'stroke-linejoin':'round', 'stroke-linecap':'round'}, svg);
      const hx = x(s.half_pg), hy = y(0.5);
      el('circle', {cx:hx, cy:hy, r:4.5, fill:cols[i % 4], stroke:P.surface, 'stroke-width':2}, svg);
      txt(svg, hx + 9, hy + 4 + (i - 1.5) * 13, s.name + ' · K_D ' + fmt(s.kd_nM, 1) + ' nM',
          {'font-size':11, fill:P.ink2});
    });

    /* crosshair: the reader aims at a concentration, never at a 2px line */
    const cross = el('line', {y1:m.t, y2:H - m.b, stroke:P.axis, 'stroke-width':1,
                              opacity:0, 'pointer-events':'none'}, svg);
    const dots = S.map((s, i) => el('circle', {r:4, fill:cols[i % 4], stroke:P.surface,
                                               'stroke-width':2, opacity:0, 'pointer-events':'none'}, svg));
    const hit = el('rect', {x:m.l, y:m.t, width:W - m.r - m.l, height:H - m.b - m.t,
                            fill:'transparent'}, svg);
    const invX = px => Math.pow(10, 1 + (px - m.l) / (W - m.r - m.l) * 6);
    function move(ev) {
      const r = svg.getBoundingClientRect(), sc = W / r.width;
      const px = (ev.clientX - r.left) * sc;
      if (px < m.l || px > W - m.r) return leave();
      const pg = invX(px);
      cross.setAttribute('x1', px); cross.setAttribute('x2', px); cross.setAttribute('opacity', 1);
      const rows = [{value: fmt(pg, 0) + ' pg/mL', name:'IL-6'}];
      S.forEach((s, i) => {
        const k = Math.max(0, Math.min(s.pts.length - 1,
          Math.round((Math.log10(pg) - 1) / 6 * (s.pts.length - 1))));
        const p = s.pts[k];
        dots[i].setAttribute('cx', x(p[0])); dots[i].setAttribute('cy', y(p[1]));
        dots[i].setAttribute('opacity', 1);
        rows.push({color:cols[i % 4], value:(p[1] * 100).toFixed(2) + '%', name:s.name + ' bound'});
      });
      showTip(ev, rows);
    }
    function leave() {
      cross.setAttribute('opacity', 0); dots.forEach(d => d.setAttribute('opacity', 0)); hideTip();
    }
    hit.addEventListener('pointermove', move);
    hit.addEventListener('pointerleave', leave);

    /* legend */
    const lg = document.createElement('div'); lg.className = 'legend';
    S.forEach((s, i) => {
      const k = document.createElement('span'); k.className = 'k';
      const ln = document.createElement('span'); ln.className = 'ln';
      ln.style.background = cols[i % 4]; k.appendChild(ln);
      k.appendChild(document.createTextNode(s.name));
      lg.appendChild(k);
    });
    const wb = document.createElement('span'); wb.className = 'k';
    const ws = document.createElement('span'); ws.className = 'sw';
    ws.style.background = P.s1; ws.style.opacity = 0.25; wb.appendChild(ws);
    wb.appendChild(document.createTextNode('wound-fluid range (post-surgical → infected)'));
    lg.appendChild(wb);
    node.appendChild(lg);
  });

  /* table-view twin */
  guard('occ-table', node => {
    const rows = DATA.occupancyTable, ps = DATA.parents || [];
    if (!rows || !rows.length || !ps.length) return empty(node, 'not yet generated');
    node.textContent = '';
    const t = document.createElement('table');
    const thead = document.createElement('thead'), tr = document.createElement('tr');
    ['IL-6 level', 'pg/mL', 'pM'].concat(ps.map(p => p.name + ' bound (%)')).forEach((h, i) => {
      const th = document.createElement('th'); th.textContent = h;
      if (!i) th.className = 'l'; tr.appendChild(th);
    });
    thead.appendChild(tr); t.appendChild(thead);
    const tb = document.createElement('tbody');
    rows.forEach(r => {
      const trr = document.createElement('tr');
      const c = [[r.label, 'l'], [fmt(r.pg_ml, 0), ''], [fmt(r.pM, 1), '']]
        .concat(ps.map(p => [(100 * (r.frac[p.name] || 0)).toFixed(2), '']));
      c.forEach(([v, cls]) => { const td = document.createElement('td');
        if (cls) td.className = cls; td.textContent = v; trr.appendChild(td); });
      tb.appendChild(trr);
    });
    t.appendChild(tb); node.appendChild(t);
  });
}

/* =========================================================================
   VIEWS 3-6 - everything scoped by the one filter row
   ========================================================================= */
const SW = DATA.switches || [];
const COLS = [
  {k:'rank', h:'Rank', d:0, l:false},
  {k:'ds', h:'DS (5′→3′)', text:true},
  {k:'window', h:'Window (0-idx)', text:true},
  {k:'ds_len', h:'DS length (nt)', d:0},
  {k:'linker_len', h:'Linker (nt)', d:0},
  {k:'gc', h:'GC (fraction)', d:2},
  {k:'dg_switch', h:'ΔG_switch (kcal/mol)', d:2},
  {k:'closed_frac', h:'Closed (fraction)', d:3},
  {k:'kd_app_nM', h:'K_D,app (nM)', d:1},
  {k:'engagement', h:'Engagement (fraction)', d:2},
  {k:'rand_covered', h:'Randomised nt covered', d:0}
];
const FILTERS = [
  {k:'dg_switch', label:'ΔG_switch (kcal/mol)', step:0.05},
  {k:'kd_app_nM', label:'K_D,app (nM)', step:0.1},
  {k:'engagement', label:'Engagement (fraction)', step:0.01}
];
const RANGE = {};
FILTERS.forEach(f => {
  const v = SW.map(r => r[f.k]);
  RANGE[f.k] = SW.length ? [Math.min.apply(null, v), Math.max.apply(null, v)] : [0, 1];
});
let sortKey = 'rank', sortDir = 1, filtered = SW.slice();

function applyFilters() {
  filtered = SW.filter(r => FILTERS.every(f => {
    const a = document.getElementById('lo-' + f.k), b = document.getElementById('hi-' + f.k);
    if (!a || !b) return true;                     /* filter row absent: show everything */
    const lo = parseFloat(a.value), hi = parseFloat(b.value);
    return (isNaN(lo) || r[f.k] >= lo) && (isNaN(hi) || r[f.k] <= hi);
  }));
  const c = document.getElementById('sw-count');
  if (c) c.textContent = filtered.length + ' of ' + SW.length + ' constructs shown';
  renderSwitchTable(); renderScatter(); renderHeatmap(); renderCoverage();
  buildStructurePicker();   /* view 7 sits below the filter row, so it follows the slice */
}
function buildFilterRow() {
  const node = card('sw-filters'); if (!node) return;
  if (!SW.length) { empty(node, 'switches.csv not yet generated — run python aptamer/switch_library.py'); return; }
  node.textContent = '';
  FILTERS.forEach(f => {
    const d = document.createElement('div'); d.className = 'f';
    const lab = document.createElement('label'); lab.textContent = f.label; d.appendChild(lab);
    const row = document.createElement('div'); row.className = 'row';
    ['lo', 'hi'].forEach((side, i) => {
      const inp = document.createElement('input');
      inp.type = 'number'; inp.step = f.step; inp.id = side + '-' + f.k;
      inp.value = RANGE[f.k][i]; inp.setAttribute('aria-label', f.label + ' ' + (i ? 'max' : 'min'));
      inp.addEventListener('change', applyFilters);
      row.appendChild(inp);
      if (!i) row.appendChild(document.createTextNode('–'));
    });
    d.appendChild(row); node.appendChild(d);
  });
  const b = document.createElement('button');
  b.textContent = 'Reset'; b.addEventListener('click', () => {
    FILTERS.forEach(f => { document.getElementById('lo-' + f.k).value = RANGE[f.k][0];
                           document.getElementById('hi-' + f.k).value = RANGE[f.k][1]; });
    sortKey = 'rank'; sortDir = 1; applyFilters();
  });
  node.appendChild(b);
  const c = document.createElement('span'); c.className = 'count'; c.id = 'sw-count';
  node.appendChild(c);
}

/* --- VIEW 3: switch library table --- */
function renderSwitchTable() {
  guard('sw-table', node => {
    if (!SW.length) return empty(node, 'switches.csv not yet generated');
    const rows = filtered.slice().sort((a, b) => {
      const x = a[sortKey], y = b[sortKey];
      if (typeof x === 'string') return sortDir * String(x).localeCompare(String(y));
      return sortDir * (x - y);
    });
    node.textContent = '';
    const wrap = document.createElement('div'); wrap.className = 'scroll';
    const t = document.createElement('table');
    const thead = document.createElement('thead'), tr = document.createElement('tr');
    COLS.forEach(c => {
      const th = document.createElement('th');
      th.className = 'sortable' + (c.text || c.k === 'rank' ? ' l' : '');
      th.textContent = c.h;
      const a = document.createElement('span'); a.className = 'arrow';
      a.textContent = sortKey === c.k ? (sortDir > 0 ? ' ▲' : ' ▼') : '';
      th.appendChild(a);
      th.title = 'sort by ' + c.h;
      th.addEventListener('click', () => {
        if (sortKey === c.k) sortDir = -sortDir;
        else { sortKey = c.k; sortDir = (c.k === 'rank' || c.text) ? 1 : 1; }
        renderSwitchTable();
      });
      tr.appendChild(th);
    });
    thead.appendChild(tr); t.appendChild(thead);
    const tb = document.createElement('tbody');
    rows.forEach(r => {
      const trr = document.createElement('tr');
      COLS.forEach(c => {
        const td = document.createElement('td');
        if (c.text || c.k === 'rank') td.className = 'l';
        if (c.text) { td.textContent = r[c.k]; td.classList.add('mono'); }
        else td.textContent = fmt(r[c.k], c.d);
        trr.appendChild(td);
      });
      trr.title = r.construct;
      tb.appendChild(trr);
    });
    t.appendChild(tb); wrap.appendChild(t); node.appendChild(wrap);
    if (!rows.length) empty(node, 'no constructs match the filters');
  });
}

/* --- VIEW 4: gain-vs-affinity scatter --- */
function renderScatter() {
  guard('sw-scatter', node => {
    if (!SW.length) return empty(node, 'switches.csv not yet generated');
    node.textContent = '';
    const W = Math.max(480, Math.min(node.clientWidth || 720, 900)), H = 400;
    const m = {t:16, r:20, b:48, l:60};
    const svg = el('svg', {width:'100%', viewBox:'0 0 ' + W + ' ' + H, role:'img',
      'aria-label':'Apparent K_D versus switch stabilisation energy'}, node);
    const dg = DATA.dgWindow;
    const x = lin(dg[0] - 0.35, dg[1] + 0.35, m.l, W - m.r);
    const ky = SW.map(r => r.kd_app_nM);
    const y0 = Math.pow(10, Math.floor(Math.log10(Math.min.apply(null, ky))));
    const y1 = Math.pow(10, Math.ceil(Math.log10(Math.max.apply(null, ky))));
    const y = log(y0, y1, H - m.b, m.t);

    /* design window as a shaded region */
    el('rect', {x:x(dg[0]), y:m.t, width:x(dg[1]) - x(dg[0]), height:H - m.b - m.t,
                fill:P.s1, opacity:0.08}, svg);
    txt(svg, (x(dg[0]) + x(dg[1])) / 2, m.t + 13,
        (DATA.dgWindowCaption || 'design window') + ' · ' +
        dg[0].toFixed(2) + ' to ' + dg[1].toFixed(2) + ' kcal/mol',
        {'text-anchor':'middle', 'font-size':10.5, fill:P.ink2});

    for (let e = Math.log10(y0); e <= Math.log10(y1) + 1e-9; e++) {
      const v = Math.pow(10, e);
      grid(svg, m.l, W - m.r, y(v), false);
      txt(svg, m.l - 8, y(v) + 4, fmt(v, 0), {'text-anchor':'end', 'font-size':10.5, fill:P.muted});
      for (let k = 2; k < 10; k++) {
        const vv = v * k; if (vv > y1) break;
        grid(svg, m.l, W - m.r, y(vv), false);
      }
    }
    for (let v = Math.ceil((dg[0] - 0.35) * 4) / 4; v <= dg[1] + 0.35; v += 0.25) {
      txt(svg, x(v), H - m.b + 16, v.toFixed(2), {'text-anchor':'middle', 'font-size':10.5, fill:P.muted});
    }
    el('line', {x1:m.l, x2:W - m.r, y1:H - m.b, y2:H - m.b, stroke:P.axis, 'shape-rendering':'crispEdges'}, svg);
    txt(svg, (m.l + W - m.r) / 2, H - 10, 'ΔG_switch (kcal/mol) — more negative = more strongly closed',
        {'text-anchor':'middle', 'font-size':11.5, fill:P.ink2});
    txt(svg, m.l - 8, m.t - 2, 'K_D,app (nM, log)', {'text-anchor':'end', 'font-size':11.5, fill:P.ink2});

    const rc = SW.map(r => r.rand_covered);
    const rMax = Math.max.apply(null, rc), rMin = Math.min.apply(null, rc);
    const eMin = Math.min.apply(null, SW.map(r => r.engagement));
    const eMax = Math.max.apply(null, SW.map(r => r.engagement));
    const rad = e => 4 + (eMax > eMin ? (e - eMin) / (eMax - eMin) : 0.5) * 4;
    const pts = [];
    filtered.forEach(r => {
      const cx = x(r.dg_switch), cy = y(r.kd_app_nM);
      const col = ordColor(r.rand_covered - rMin, rMax - rMin + 1);
      el('circle', {cx:cx, cy:cy, r:rad(r.engagement), fill:col, 'fill-opacity':0.85,
                    stroke:P.surface, 'stroke-width':2}, svg);
      pts.push({x:cx, y:cy, r:r, col:col});
    });

    /* nearest-point hover: the pointer only has to be closest, not dead centre */
    const halo = el('circle', {r:0, fill:'none', stroke:P.ink, 'stroke-width':1.5, opacity:0}, svg);
    const hit = el('rect', {x:m.l, y:m.t, width:W - m.r - m.l, height:H - m.b - m.t, fill:'transparent'}, svg);
    hit.addEventListener('pointermove', ev => {
      const bb = svg.getBoundingClientRect(), sc = W / bb.width;
      const px = (ev.clientX - bb.left) * sc, py = (ev.clientY - bb.top) * sc;
      let best = null, bd = 1e9;
      pts.forEach(p => { const d = (p.x - px) ** 2 + (p.y - py) ** 2;
                         if (d < bd) { bd = d; best = p; } });
      if (!best || bd > 900) { halo.setAttribute('opacity', 0); return hideTip(); }
      halo.setAttribute('cx', best.x); halo.setAttribute('cy', best.y);
      halo.setAttribute('r', rad(best.r.engagement) + 4); halo.setAttribute('opacity', 0.9);
      const r = best.r;
      showTip(ev, [
        {value:'#' + r.rank + ' ' + r.ds, color:best.col},
        {value:fmt(r.kd_app_nM, 1) + ' nM', name:'K_D,app'},
        {value:r.dg_switch.toFixed(2) + ' kcal/mol', name:'ΔG_switch'},
        {value:r.engagement.toFixed(2), name:'engagement'},
        {value:String(r.rand_covered), name:'randomised nt covered'},
        {value:r.window, name:'window · linker ' + r.linker_len + ' nt'}
      ]);
    });
    hit.addEventListener('pointerleave', () => { halo.setAttribute('opacity', 0); hideTip(); });

    /* legends: colour ramp + size key */
    const lg = document.createElement('div'); lg.className = 'legend';
    const ramp = document.createElement('span'); ramp.className = 'rampbar';
    ramp.appendChild(document.createTextNode('randomised nt covered ' + rMin));
    const bar = document.createElement('span'); bar.className = 'bar';
    const stops = [];
    for (let i = rMin; i <= rMax; i++) stops.push(ordColor(i - rMin, rMax - rMin + 1));
    bar.style.background = 'linear-gradient(90deg,' + stops.join(',') + ')';
    ramp.appendChild(bar);
    ramp.appendChild(document.createTextNode(String(rMax)));
    lg.appendChild(ramp);
    const sz = document.createElement('span'); sz.className = 'k';
    sz.textContent = 'point size = engagement (' + eMin.toFixed(2) + '–' + eMax.toFixed(2) + ')';
    lg.appendChild(sz);
    node.appendChild(lg);
  });
}

/* --- VIEW 5: linker x DS-length heatmap --- */
function renderHeatmap() {
  guard('sw-heatmap', node => {
    if (!SW.length) return empty(node, 'switches.csv not yet generated');
    node.textContent = '';
    const dsL = [...new Set(SW.map(r => r.ds_len))].sort((a, b) => a - b);
    const lnL = [...new Set(SW.map(r => r.linker_len))].sort((a, b) => a - b);
    const best = {}, cnt = {};
    filtered.forEach(r => {
      const k = r.ds_len + '|' + r.linker_len;
      if (best[k] === undefined || r.kd_app_nM < best[k]) best[k] = r.kd_app_nM;
      cnt[k] = (cnt[k] || 0) + 1;
    });
    /* fixed global domain: filtering must not repaint the survivors */
    const all = SW.map(r => r.kd_app_nM);
    const lo = Math.min.apply(null, all), hi = Math.max.apply(null, all);
    const cw = 52, ch = 34, gap = 2, m = {t:24, r:16, b:34, l:96};
    const W = m.l + lnL.length * cw + m.r, H = m.t + dsL.length * ch + m.b;
    const svg = el('svg', {width:'100%', viewBox:'0 0 ' + W + ' ' + H, role:'img',
      'aria-label':'Best apparent K_D for each linker length and DS length'}, node);
    lnL.forEach((ln, j) => txt(svg, m.l + j * cw + cw / 2, m.t - 8, String(ln),
      {'text-anchor':'middle', 'font-size':10.5, fill:P.muted}));
    txt(svg, m.l + lnL.length * cw / 2, H - 8, 'linker length (nt)',
        {'text-anchor':'middle', 'font-size':11.5, fill:P.ink2});
    txt(svg, m.l - 10, m.t - 8, 'DS length (nt)', {'text-anchor':'end', 'font-size':11.5, fill:P.ink2});
    dsL.forEach((ds, i) => {
      txt(svg, m.l - 10, m.t + i * ch + ch / 2 + 4, String(ds),
          {'text-anchor':'end', 'font-size':10.5, fill:P.muted});
      lnL.forEach((ln, j) => {
        const k = ds + '|' + ln, v = best[k];
        const gx = m.l + j * cw + gap / 2, gy = m.t + i * ch + gap / 2;
        if (v === undefined) {
          el('rect', {x:gx, y:gy, width:cw - gap, height:ch - gap, rx:4, fill:'none',
                      stroke:P.grid, 'stroke-width':1}, svg);
          return;
        }
        const t = hi > lo ? (v - lo) / (hi - lo) : 0.5;
        const col = seqColor(t);
        const cell = el('rect', {x:gx, y:gy, width:cw - gap, height:ch - gap, rx:4, fill:col}, svg);
        txt(svg, gx + (cw - gap) / 2, gy + (ch - gap) / 2 + 4, fmt(v, 0),
            {'text-anchor':'middle', 'font-size':10.5, fill:inkOn(col)});
        cell.addEventListener('pointermove', ev => {
          cell.setAttribute('stroke', P.ink); cell.setAttribute('stroke-width', 1.5);
          showTip(ev, [{value:fmt(v, 1) + ' nM', name:'best K_D,app', color:col},
                       {value:'DS ' + ds + ' nt · linker ' + ln + ' nt'},
                       {value:cnt[k] + ' construct(s)', name:'in this cell'}]);
        });
        cell.addEventListener('pointerleave', () => { cell.removeAttribute('stroke'); hideTip(); });
      });
    });
    const lg = document.createElement('div'); lg.className = 'legend';
    const ramp = document.createElement('span'); ramp.className = 'rampbar';
    ramp.appendChild(document.createTextNode('best K_D,app  ' + fmt(lo, 0) + ' nM'));
    const bar = document.createElement('span'); bar.className = 'bar';
    bar.style.background = 'linear-gradient(90deg,' + P.seq.join(',') + ')';
    ramp.appendChild(bar);
    ramp.appendChild(document.createTextNode(fmt(hi, 0) + ' nM'));
    lg.appendChild(ramp);
    const nk = document.createElement('span'); nk.className = 'k';
    const sw = document.createElement('span'); sw.className = 'sw';
    sw.style.background = 'transparent'; sw.style.border = '1px solid ' + P.grid;
    nk.appendChild(sw); nk.appendChild(document.createTextNode('no passing construct'));
    lg.appendChild(nk);
    node.appendChild(lg);
  });
}

/* --- VIEW 6: window coverage track --- */
function renderCoverage() {
  guard('sw-coverage', node => {
    const seq = DATA.parentSeq;
    if (!SW.length || !seq) return empty(node, 'needs both parents.json and switches.csv');
    node.textContent = '';
    const n = seq.length, rand = new Set(DATA.randomised);
    const cov = new Array(n).fill(0);
    filtered.forEach(r => { for (let i = r.w0; i <= r.w1 && i < n; i++) cov[i]++; });
    const maxc = Math.max(1, Math.max.apply(null, cov));
    const cellW = 15, m = {t:18, r:16, b:64, l:44};
    const barH = 150;
    const W = m.l + n * cellW + m.r, H = m.t + barH + m.b;
    const svg = el('svg', {width:'100%', viewBox:'0 0 ' + W + ' ' + H, role:'img',
      'aria-label':'Constructs covering each position of the parent aptamer'}, svg0(node));
    const y = lin(0, maxc, m.t + barH, m.t);
    [0, 0.5, 1].forEach(f => {
      const v = maxc * f;
      grid(svg, m.l, W - m.r, y(v), f === 0);
      txt(svg, m.l - 8, y(v) + 4, fmt(v, 0), {'text-anchor':'end', 'font-size':10.5, fill:P.muted});
    });
    txt(svg, m.l - 8, m.t - 4, 'constructs', {'text-anchor':'end', 'font-size':11, fill:P.ink2});

    /* randomised library blocks: the wash marks them behind everything */
    let i = 0;
    while (i < n) {
      if (rand.has(i)) {
        let j = i; while (j + 1 < n && rand.has(j + 1)) j++;
        el('rect', {x:m.l + i * cellW, y:m.t, width:(j - i + 1) * cellW, height:barH + 34,
                    fill:P.s2, opacity:0.13, rx:3}, svg);
        i = j + 1;
      } else i++;
    }
    /* coverage columns */
    for (let k = 0; k < n; k++) {
      const bw = cellW - 3;
      const gx = m.l + k * cellW + 1.5;
      const h = y(0) - y(cov[k]);
      if (h > 0.5) el('rect', {x:gx, y:y(cov[k]), width:bw, height:h, rx:Math.min(4, bw / 2),
                               fill:P.s1}, svg);
      const t = txt(svg, gx + bw / 2, m.t + barH + 18, seq[k],
                    {'text-anchor':'middle', 'font-size':10.5,
                     fill:rand.has(k) ? P.ink : P.muted,
                     'font-family':'ui-monospace, Consolas, monospace'});
      if (rand.has(k)) t.setAttribute('font-weight', '600');
      if (k % 10 === 0) txt(svg, gx + bw / 2, m.t + barH + 32, String(k),
                            {'text-anchor':'middle', 'font-size':9.5, fill:P.muted});
      const hit = el('rect', {x:m.l + k * cellW, y:m.t, width:cellW, height:barH + 24,
                              fill:'transparent'}, svg);
      hit.addEventListener('pointermove', ev => showTip(ev, [
        {value:cov[k] + ' construct(s)', name:'target position ' + k, color:P.s1},
        {value:seq[k] + ' · ' + (rand.has(k) ? 'randomised block' : 'fixed region')}
      ]));
      hit.addEventListener('pointerleave', hideTip);
    }
    txt(svg, m.l + n * cellW / 2, H - 22, DATA.switchParent + ' — position in the ' + n + '-nt parent (0-indexed)',
        {'text-anchor':'middle', 'font-size':11.5, fill:P.ink2});

    const lg = document.createElement('div'); lg.className = 'legend';
    [[P.s1, 'constructs whose DS window covers this position', 'sw'],
     ['rgba(235,104,52,0.35)', 'randomised library block (bold letters)', 'sw']]
      .forEach(([c, label]) => {
        const k = document.createElement('span'); k.className = 'k';
        const s = document.createElement('span'); s.className = 'sw'; s.style.background = c;
        k.appendChild(s); k.appendChild(document.createTextNode(label)); lg.appendChild(k);
      });
    node.appendChild(lg);
  });
}
/* the coverage track is wider than the card; give it its own horizontal scroller */
function svg0(node) {
  const d = document.createElement('div');
  d.style.overflowX = 'auto'; node.appendChild(d); return d;
}

/* =========================================================================
   VIEW 7 - secondary structure: fornac 2D layout, arc diagram as the fallback
   ========================================================================= */
function pairs(db) {
  const st = [], out = [];
  for (let i = 0; i < db.length; i++) {
    const c = db[i];
    if (c === '(') st.push(i);
    else if (c === ')') { const j = st.pop(); if (j !== undefined) out.push([j, i]); }
  }
  return out;
}
function structureOptions() {
  const opts = [];
  (DATA.parents || []).forEach(p => {
    opts.push({label:'parent ' + p.name + ' @ 22 °C (MFE ' + p.mfe_22C.toFixed(2) + ')',
               seq:p.sequence, db:p.structure_22C, kind:'parent'});
    opts.push({label:'parent ' + p.name + ' @ 37 °C (MFE ' + p.mfe_37C.toFixed(2) + ')',
               seq:p.sequence, db:p.structure_37C, kind:'parent'});
  });
  /* switch constructs follow the filter row, like every other view below it */
  filtered.slice(0, 25).forEach(r => opts.push({
    label:'switch #' + r.rank + ' · DS ' + r.ds + ' · window ' + r.window +
          ' · linker ' + r.linker_len,
    seq:r.construct, db:r.structure, kind:'switch', row:r
  }));
  return opts;
}
/* ---------- fornac (primary renderer) ----------
   fornac 1.1.8 is a webpack UMD bundle, so in a browser it lands on window.fornac; the
   forna site's own build also exposes window.FornaContainer. Accept either, and treat
   "neither is there" as the signal to fall back. d3 3.5.13 is bundled inside fornac.js
   and jQuery is never referenced, so no other CDN script is needed. */
function fornaCtor() {
  return window.FornaContainer ||
         (window.fornac && window.fornac.FornaContainer) || null;
}
/* fornac carries its own stylesheet and style-loader injects it into <head> unscoped. Its
   bare `svg { width:100%; min-height:100% }` and `text {}` rules would resize every other
   chart on this page, so delete that copy -- a copy scoped to #ss-forna is already inlined
   in the page's own <style data-dash>. */
function dropFornacGlobalCss() {
  const styles = document.querySelectorAll('head style:not([data-dash])');
  for (let i = 0; i < styles.length; i++) {
    const t = styles[i].textContent || '';
    if (t.indexOf('circle.outline_node') >= 0 && /(^|\})\s*svg\s*\{/.test(t))
      styles[i].remove();
  }
}
/* [start, end) regions of the construct, in the page palette */
function ssRegions(opt) {
  const n = opt.seq.length;
  if (opt.kind !== 'switch') return [[0, n, 'aptamer', P.s1]];
  const nDs = opt.row.ds_len, nLink = opt.row.linker_len, nAp = n - nDs - nLink;
  return [[0, nAp, 'aptamer', P.s1], [nAp, nAp + nLink, 'poly-T linker', P.s3],
          [nAp + nLink, n, 'displacement strand', P.s2]];
}
function ssLegend(node, keys, shape) {
  const lg = document.createElement('div'); lg.className = 'legend';
  keys.forEach(([c, label]) => {
    const k = document.createElement('span'); k.className = 'k';
    if (c) { const s = document.createElement('span'); s.className = shape || 'sw';
             s.style.background = c; k.appendChild(s); }
    k.appendChild(document.createTextNode(label)); lg.appendChild(k);
  });
  node.appendChild(lg);
}
function fornaView(node, opt) {
  const Forna = fornaCtor(), seq = opt.seq, db = opt.db, n = seq.length;
  const host = document.createElement('div');
  host.id = 'ss-forna';                       /* the id every fornac CSS rule is scoped to */
  const W = Math.max(320, node.clientWidth || 720);
  const H = Math.max(360, Math.min(560, 200 + 3 * n));
  host.style.height = H + 'px';
  node.appendChild(host);
  dropFornacGlobalCss();
  const c = new Forna('#ss-forna', {initialSize:[W, H], labelInterval:10,
    allowPanningAndZooming:true, applyForce:true, transitionDuration:0});
  /* colour every nucleotide by construct region. fornac's customColors is
     {colorValues: {structureName: {1-based position: colour}}}, and '' matches any
     structure, so the RNA does not have to be named. */
  const cv = {}, regions = ssRegions(opt);
  regions.forEach(([a, b, label, col]) => { for (let i = a; i < b; i++) cv[i + 1] = col; });
  c.addCustomColors({colorValues: {'': cv}});
  /* the sequence is DNA and is passed through verbatim -- no T->U conversion, so the
     letters on screen are the letters you would order. fornac only needs the dot-bracket
     for the layout and never validates the alphabet. */
  c.addRNA(db, {sequence: seq, labelInterval: 10});
  c.changeColorScheme('custom');
  ssLegend(node, regions.map(([a, b, label, col]) => [col, label]).concat(
    [[null, pairs(db).length + ' base pairs · ' + n + ' nt'],
     [null, 'drag to pan · scroll to zoom · drag a nucleotide to pull the layout apart']]));
}
/* ---------- base-pair arc diagram (offline fallback) ---------- */
function arcView(node, opt, reason) {
  if (reason) {
    const w = document.createElement('div'); w.className = 'issues';
    w.style.margin = '0 0 12px'; w.textContent = reason; node.appendChild(w);
  }
  const seq = opt.seq, db = opt.db, n = seq.length, pr = pairs(db);
  const cellW = 14, m = {t:118, l:26, r:26, b:60};
  const W = m.l + n * cellW + m.r, H = m.t + m.b;
  const host = document.createElement('div'); host.style.overflowX = 'auto';
  node.appendChild(host);
  const svg = el('svg', {width:'100%', viewBox:'0 0 ' + W + ' ' + H, role:'img',
    'aria-label':'Base-pair arc diagram'}, host);
  const base = m.t, X = i => m.l + i * cellW + cellW / 2;

  /* construct regions: aptamer | poly-T linker | displacement strand */
  let nAp = n, nLink = 0, nDs = 0;
  if (opt.kind === 'switch') { nDs = opt.row.ds_len; nLink = opt.row.linker_len; nAp = n - nDs - nLink; }
  const regions = opt.kind === 'switch'
    ? [[0, nAp, 'aptamer', P.s1], [nAp, nAp + nLink, 'poly-T linker', P.s3],
       [nAp + nLink, n, 'displacement strand', P.s2]]
    : [[0, n, 'aptamer', P.s1]];
  regions.forEach(([a, b, label, col]) => {
    el('rect', {x:X(a) - cellW / 2 + 1, y:base + 20, width:(b - a) * cellW - 2, height:5,
                rx:2.5, fill:col}, svg);
    txt(svg, (X(a) + X(b - 1)) / 2, base + 40, label,
        {'text-anchor':'middle', 'font-size':10.5, fill:P.ink2});
  });
  /* randomised blocks on the parent portion */
  const rand = new Set(DATA.randomised);
  /* arcs above the sequence line */
  const maxSpan = pr.reduce((s, p) => Math.max(s, p[1] - p[0]), 1);
  pr.forEach(([a, b]) => {
    const x1 = X(a), x2 = X(b), r = (b - a) / maxSpan;
    const h = 18 + r * (m.t - 34);
    const cross = opt.kind === 'switch' && ((a < nAp) !== (b < nAp));
    el('path', {d:'M' + x1 + ' ' + base + ' C ' + x1 + ' ' + (base - h) + ' ' +
                   x2 + ' ' + (base - h) + ' ' + x2 + ' ' + base,
                fill:'none', stroke:cross ? P.s2 : P.s1, 'stroke-width':2,
                opacity:0.85, 'stroke-linecap':'round'}, svg);
  });
  el('line', {x1:m.l, x2:W - m.r, y1:base, y2:base, stroke:P.axis, 'shape-rendering':'crispEdges'}, svg);
  for (let i = 0; i < n; i++) {
    const t = txt(svg, X(i), base + 14, seq[i], {'text-anchor':'middle', 'font-size':10.5,
      fill:(i < nAp && rand.has(i)) ? P.ink : P.muted,
      'font-family':'ui-monospace, Consolas, monospace'});
    if (i < nAp && rand.has(i)) t.setAttribute('font-weight', '600');
  }
  /* per-position hover */
  const partner = new Array(n).fill(-1);
  pr.forEach(([a, b]) => { partner[a] = b; partner[b] = a; });
  for (let i = 0; i < n; i++) {
    const hit = el('rect', {x:X(i) - cellW / 2, y:8, width:cellW, height:base + 10,
                            fill:'transparent'}, svg);
    hit.addEventListener('pointermove', ev => showTip(ev, [
      {value:seq[i] + ' at ' + i, color:P.s1},
      {value:partner[i] >= 0 ? 'paired with ' + seq[partner[i]] + ' at ' + partner[i] : 'unpaired'},
      {value:(i < nAp ? (rand.has(i) ? 'aptamer, randomised block' : 'aptamer, fixed region')
                      : (i < nAp + nLink ? 'poly-T linker' : 'displacement strand'))}
    ]));
    hit.addEventListener('pointerleave', hideTip);
  }

  const keys = opt.kind === 'switch'
    ? [[P.s1, 'base pair within the aptamer'], [P.s2, 'DS ↔ aptamer pair (the switch)']]
    : [[P.s1, 'base pair']];
  ssLegend(node, keys.concat([[null, pr.length + ' pairs · ' + n + ' nt']]), 'ln');
}
/* the dot-bracket itself, under whichever renderer drew above it */
function textView(node, opt) {
  const det = document.createElement('details'); det.className = 'tv';
  const sm = document.createElement('summary');
  sm.textContent = 'dot-bracket + sequence (text view)';
  det.appendChild(sm);
  const pre = document.createElement('div'); pre.className = 'seq';
  pre.style.marginTop = '6px';
  pre.textContent = opt.seq + '\n' + opt.db;
  pre.style.whiteSpace = 'pre-wrap';
  det.appendChild(pre); node.appendChild(det);
}
const FORNA_MISSING =
  'fornac is not available — either its CDN script did not load or the generator could ' +
  'not read aptamer/vendor/fornac.css — so this panel is showing the built-in base-pair ' +
  'arc diagram instead. The structure is identical; only the layout differs. Every other ' +
  'view except the 3D panel works offline either way.';
function renderStructure(opt) {
  guard('ss-chart', node => {
    node.textContent = '';
    if (!opt) return empty(node, 'no structures available yet');
    if (fornaCtor()) {
      try {
        fornaView(node, opt);
      } catch (e) {
        node.textContent = '';
        arcView(node, opt, 'fornac loaded but could not draw this structure (' + e.message +
                           '); showing the base-pair arc diagram instead.');
      }
    } else {
      arcView(node, opt, FORNA_MISSING);
    }
    textView(node, opt);
  });
}
let ssOpts = [];
function buildStructurePicker() {
  const sel = document.getElementById('ss-select'); if (!sel) return;
  const keep = sel.selectedIndex >= 0 ? sel.options[sel.selectedIndex].textContent : null;
  ssOpts = structureOptions();
  if (!ssOpts.length) { sel.textContent = ''; renderStructure(null); return; }
  sel.textContent = '';
  ssOpts.forEach((o, i) => { const op = document.createElement('option');
    op.value = i; op.textContent = o.label; sel.appendChild(op); });
  const again = ssOpts.findIndex(o => o.label === keep);
  sel.value = again >= 0 ? again : 0;
  if (!sel.dataset.wired) {
    sel.addEventListener('change', () => renderStructure(ssOpts[+sel.value]));
    sel.dataset.wired = '1';
  }
  window.__ssRedraw = () => renderStructure(ssOpts[+sel.value || 0]);
  renderStructure(ssOpts[+sel.value || 0]);
}

/* =========================================================================
   VIEW 8 - 3D viewer (3Dmol.js, needs internet)
   ========================================================================= */
function render3D() {
  const host = document.getElementById('mol-viewer');
  const sel = document.getElementById('mol-select');
  if (!host || !sel) return;
  const structs = DATA.structures || [];
  if (!structs.length) { empty(host, 'no .pdb / .cif in aptamer/ yet — not yet generated'); return; }
  sel.textContent = '';
  structs.forEach((s, i) => { const o = document.createElement('option');
    o.value = i; o.textContent = s.name + ' — ' + s.label; sel.appendChild(o); });
  function draw() {
    if (!window.$3Dmol) {
      empty(host, '3D viewer could not load. 3Dmol.js is fetched from a CDN, so this one panel needs internet; every other view works offline.');
      return;
    }
    const s = structs[+sel.value || 0];
    host.textContent = '';
    const v = $3Dmol.createViewer(host, {backgroundColor: P.surface});
    v.addModel(s.text, s.format);
    v.setStyle({}, {cartoon:{color:P.s1, opacity:0.95}});
    /* co-folds always write the protein as chain A and the DNA as chain B, so colour by
       chain first; the residue-name rule below then catches any file that does not. */
    if (s.chainA) v.setStyle({chain:'A'}, {cartoon:{color:P.s1, opacity:0.95}});
    if (s.chainB) v.setStyle({chain:'B'}, {stick:{color:P.s2, radius:0.22},
                                           cartoon:{color:P.s2, style:'trace'}});
    v.setStyle({resn:['DA','DC','DG','DT','A','C','G','U','T']},
               {stick:{colorscheme:'default', color:P.s2, radius:0.22},
                cartoon:{color:P.s2, style:'trace'}});
    v.setStyle({hetflag:true}, {stick:{color:P.s4, radius:0.18}});
    v.zoomTo(); v.render();
    /* say which colour is which for the structure actually on screen */
    const cap = document.getElementById('mol-key');
    if (cap) {
      cap.textContent = '';
      const keys = s.chainA
        ? [[P.s1, 'chain A — ' + s.chainA], [P.s2, 'chain B — ' + s.chainB]]
        : [[P.s1, 'protein backbone'], [P.s2, 'nucleic acid, if any'],
           [P.s4, 'ligands / heteroatoms']];
      keys.forEach(([c, label]) => {
        const k = document.createElement('span'); k.className = 'k';
        const sw = document.createElement('span'); sw.className = 'sw'; sw.style.background = c;
        k.appendChild(sw); k.appendChild(document.createTextNode(label)); cap.appendChild(k);
      });
      const k = document.createElement('span'); k.className = 'k';
      k.appendChild(document.createTextNode(s.atoms.toLocaleString() + ' atoms'));
      cap.appendChild(k);
    }
    window.__molRedraw = draw;
  }
  sel.addEventListener('change', draw);
  if (window.$3Dmol) draw();
  else {
    const t = setInterval(() => { if (window.$3Dmol) { clearInterval(t); draw(); } }, 200);
    setTimeout(() => { clearInterval(t); if (!window.$3Dmol) draw(); }, 6000);
  }
}

/* =========================================================================
   VIEW 9 - mismatch refinement: before -> after at a FIXED tether
   Form: dumbbell. The job is "before -> after per item" (one hue, two shades), and the item
   identity is the variant, so a bar chart of the after-value alone would throw away the whole
   point -- that the tether did not move. Direction of change carries polarity, so the
   connector takes the diverging pair (cool = tighter, warm = looser) while the two endpoint
   marks stay one hue in two shades and are told apart by fill, not only colour.
   ========================================================================= */
const MM = DATA.mismatches || [];
const MMH = DATA.mismatchHeadline;
const MM_EMPTY = 'mismatches.csv not yet generated — run python aptamer/mismatch_tune.py';

function renderMismatchChart() {
  guard('mm-chart', node => {
    if (!MM.length) return empty(node, MM_EMPTY);
    node.textContent = '';
    const rows = MM.slice().sort((a, b) => a.d_kd_nM - b.d_kd_nM);
    const rowH = 19, m = {t:52, r:172, b:58, l:120};
    const W = Math.max(620, Math.min(node.clientWidth || 900, 1000));
    const H = m.t + rows.length * rowH + m.b;
    const svg = el('svg', {width:'100%', viewBox:'0 0 ' + W + ' ' + H, role:'img',
      'aria-label':'Apparent K_D of each construct before and after a single DS mismatch, ' +
                   'at unchanged tether length'}, node);
    const ref = MMH && MMH.matchAlt;      /* the shortest-tether unmutated alternative */
    const vals = [];
    rows.forEach(r => vals.push(r.kd_app_nM, r.kd_app_wt_nM));
    if (ref) vals.push(ref.kd_app_nM);
    const lo = Math.pow(10, Math.floor(Math.log10(Math.min.apply(null, vals))));
    const hi = Math.pow(10, Math.ceil(Math.log10(Math.max.apply(null, vals))));
    const x = log(lo, hi, m.l, W - m.r);
    const vline = (v, col, dash) => el('line', {x1:x(v), x2:x(v), y1:m.t - 6,
      y2:m.t + rows.length * rowH + 4, stroke:col, 'stroke-width':1,
      'stroke-dasharray':dash || null, 'shape-rendering':dash ? null : 'crispEdges'}, svg);
    for (let e = Math.log10(lo); e <= Math.log10(hi) + 1e-9; e++) {
      const v = Math.pow(10, e);
      vline(v, P.grid);
      txt(svg, x(v), m.t + rows.length * rowH + 20, fmt(v, 0),
          {'text-anchor':'middle', 'font-size':10.5, fill:P.muted});
      for (let k = 2; k < 10; k++) { const vv = v * k; if (vv > hi) break; vline(vv, P.grid); }
    }
    txt(svg, (m.l + W - m.r) / 2, H - 20,
        'K_D,app (nM, log) — open mark = unmutated construct, filled mark = single mismatch',
        {'text-anchor':'middle', 'font-size':11.5, fill:P.ink2});

    /* the construct the mismatch is arguing against: the shortest-tether UNMUTATED construct
       that matches the best variant's affinity while covering at least as much, and what that
       costs in tether. */
    if (ref) {
      vline(ref.kd_app_nM, P.ink2, '4 3');
      /* the annotation sits on whichever side of the rule has room */
      const mid = (m.l + W - m.r) / 2, left = x(ref.kd_app_nM) < mid;
      const ax = x(ref.kd_app_nM) + (left ? 7 : -7), anch = left ? 'start' : 'end';
      txt(svg, ax, m.t - 24,
          'no mismatch, covering ≥ ' + MMH.covered + ' randomised nt: ' +
          fmt(ref.kd_app_nM, 0) + ' nM',
          {'text-anchor':anch, 'font-size':10.5, fill:P.ink2});
      txt(svg, ax, m.t - 12, 'but it costs a ' + ref.tether_nt + '-nt tether',
          {'text-anchor':anch, 'font-size':10.5, fill:P.ink2});
    }
    /* the two right-hand columns: the tether never moves, which is the entire claim */
    txt(svg, W - m.r + 12, m.t - 24, 'tether', {'font-size':10.5, fill:P.ink2});
    txt(svg, W - m.r + 12, m.t - 12, '(nt, fixed)', {'font-size':10, fill:P.muted});
    txt(svg, W - m.r + 84, m.t - 24, 'same K_D', {'font-size':10.5, fill:P.ink2});
    txt(svg, W - m.r + 84, m.t - 12, 'via linker (nt)', {'font-size':10, fill:P.muted});

    const best3 = rows.slice().sort((a, b) => a.kd_app_nM - b.kd_app_nM).slice(0, 3)
                      .map(r => r.name);
    rows.forEach((r, i) => {
      const y = m.t + i * rowH + rowH / 2;
      const better = r.d_kd_nM < 0;
      const col = better ? P.s1 : P.s2;
      el('line', {x1:x(r.kd_app_wt_nM), x2:x(r.kd_app_nM), y1:y, y2:y, stroke:col,
                  'stroke-width':2, 'stroke-linecap':'round'}, svg);
      el('circle', {cx:x(r.kd_app_wt_nM), cy:y, r:4.5, fill:P.surface, stroke:P.muted,
                    'stroke-width':2}, svg);
      el('circle', {cx:x(r.kd_app_nM), cy:y, r:5, fill:col, stroke:P.surface,
                    'stroke-width':2}, svg);
      txt(svg, m.l - 10, y + 3.5, r.name, {'text-anchor':'end', 'font-size':10.5,
        fill:P.ink2, 'font-family':'ui-monospace, Consolas, monospace'});
      txt(svg, W - m.r + 12, y + 3.5, String(r.tether_nt),
          {'font-size':10.5, fill:P.ink2});
      txt(svg, W - m.r + 84, y + 3.5,
          r.tether_equiv === null || r.tether_equiv === undefined ? '—' : String(r.tether_equiv),
          {'font-size':10.5, fill:P.muted});
      /* direct-label only the three tightest, not every mark */
      if (best3.indexOf(r.name) >= 0)
        txt(svg, x(r.kd_app_nM) - 10, y + 3.5, fmt(r.kd_app_nM, 0) + ' nM',
            {'text-anchor':'end', 'font-size':10.5, fill:P.ink});
      const band = el('rect', {x:m.l, y:y - rowH / 2, width:W - m.r - m.l, height:rowH,
                               fill:'transparent'}, svg);
      band.addEventListener('pointermove', ev => showTip(ev, [
        {value:r.name, color:col},
        {value:r.ds_wt + ' → ' + r.ds, name:'DS, position ' + r.mismatch_pos},
        {value:fmt(r.kd_app_wt_nM, 1) + ' → ' + fmt(r.kd_app_nM, 1) + ' nM', name:'K_D,app'},
        {value:(r.d_kd_nM > 0 ? '+' : '') + fmt(r.d_kd_nM, 1) + ' nM',
         name:better ? 'tighter' : 'looser'},
        {value:r.tether_nt + ' nt (unchanged)', name:'tether'},
        {value:r.tether_equiv === null || r.tether_equiv === undefined ? 'not reachable' :
               r.tether_equiv + ' nt', name:'tether needed to buy this K_D by lengthening'},
        {value:r.dg_switch.toFixed(2) + ' kcal/mol', name:'ΔG_switch'},
        {value:r.closed_frac.toFixed(3), name:'closed fraction'},
        {value:r.engagement.toFixed(2), name:'engagement'},
        {value:String(r.rand_covered), name:'randomised nt covered'}
      ]));
      band.addEventListener('pointerleave', hideTip);
    });

    const lg = document.createElement('div'); lg.className = 'legend';
    const key = (label, make) => {
      const k = document.createElement('span'); k.className = 'k';
      k.appendChild(make()); k.appendChild(document.createTextNode(label)); lg.appendChild(k);
    };
    key('unmutated construct (K_D,app before)', () => {
      const s = document.createElement('span'); s.className = 'sw';
      s.style.borderRadius = '50%'; s.style.background = 'transparent';
      s.style.border = '2px solid ' + P.muted; return s;
    });
    key('single mismatch tightens K_D,app', () => {
      const s = document.createElement('span'); s.className = 'sw';
      s.style.borderRadius = '50%'; s.style.background = P.s1; return s;
    });
    key('single mismatch loosens K_D,app', () => {
      const s = document.createElement('span'); s.className = 'sw';
      s.style.borderRadius = '50%'; s.style.background = P.s2; return s;
    });
    if (ref) key('shortest-tether unmutated construct of equal affinity (dashed rule)', () => {
      const s = document.createElement('span'); s.className = 'ln';
      s.style.background = 'transparent';
      s.style.borderTop = '2px dashed ' + P.ink2; s.style.height = '0'; return s;
    });
    node.appendChild(lg);
  });
}

const MM_COLS = [
  {k:'name', h:'Variant (Wilson convention)', text:true},
  {k:'ds_wt', h:'DS unmutated (5′→3′)', text:true},
  {k:'ds', h:'DS mismatched (5′→3′)', text:true},
  {k:'mismatch_pos', h:'Mismatch position (1-idx)', d:0},
  {k:'window', h:'Window (0-idx)', text:true},
  {k:'linker_len', h:'Linker (nt)', d:0},
  {k:'tether_nt', h:'Tether (nt)', d:0, tip:'linker + DS; a mismatch never changes it'},
  {k:'dg_switch', h:'ΔG_switch (kcal/mol)', d:2},
  {k:'closed_frac', h:'Closed (fraction)', d:3},
  {k:'kd_app_wt_nM', h:'K_D,app unmutated (nM)', d:1},
  {k:'kd_app_nM', h:'K_D,app mismatched (nM)', d:1},
  {k:'d_kd_nM', h:'ΔK_D,app (nM)', d:1, cell:(td, r) => {
    const dot = document.createElement('span');
    dot.className = 'dot'; dot.style.background = r.d_kd_nM < 0 ? P.s1 : P.s2;
    td.appendChild(dot);
    td.appendChild(document.createTextNode((r.d_kd_nM > 0 ? '+' : '') + fmt(r.d_kd_nM, 1)));
    td.title = r.d_kd_nM < 0 ? 'tighter than the unmutated construct' : 'looser';
  }},
  {k:'engagement', h:'Engagement (fraction)', d:2},
  {k:'selectivity', h:'Selectivity (kcal/mol)', d:2},
  {k:'rand_covered', h:'Randomised nt covered', d:0},
  {k:'tether_equiv', h:'Tether for same K_D w/o mismatch (nt)', d:0,
   tip:'shortest tether at which an unmutated construct covering at least as many randomised ' +
       'nt reaches this K_D,app'}
];
const mmSort = {k:'d_kd_nM', d:1, draw:() => renderMismatchTable()};
function renderMismatchTable() {
  guard('mm-table', node => {
    if (!MM.length) return empty(node, MM_EMPTY);
    sortableTable(node, MM_COLS, MM, mmSort);
  });
}

/* =========================================================================
   VIEW 10 - negative controls
   Form: small-multiple histograms of engagement, one facet per arm, shared x. Faceting
   carries arm identity, which frees colour to do the job that matters -- emphasis: the
   designed arm in the accent hue, the three controls in the de-emphasis grey. A 4-series
   overlay would have made the reader hunt for the gap that is the entire result.
   ========================================================================= */
const NC = DATA.controls || [];
const NCSUM = DATA.controlSummary || [];
const NC_EMPTY =
  'negative_controls.csv not yet generated — run python aptamer/negative_controls.py';
const ENG_MIN = (DATA.budget && DATA.budget.MIN_ENGAGEMENT) || 0.6;

function renderControlDist() {
  guard('nc-chart', node => {
    if (!NC.length || !NCSUM.length) return empty(node, NC_EMPTY);
    node.textContent = '';
    const arms = NCSUM.map(s => s.arm);
    const NB = 20, bw = 1 / NB;                      /* engagement is a fraction on [0,1] */
    const counts = {};
    arms.forEach(a => { counts[a] = new Array(NB).fill(0); });
    NC.forEach(r => {
      if (!counts[r.arm]) return;
      /* `*NB` with an epsilon, not `/bw`: 0.6/0.05 is 11.999999999999998 in binary floating
         point, which would drop every construct sitting exactly ON the engagement filter into
         the bin to the LEFT of the threshold rule -- the one place a reader must not be
         misled by a rounding artefact. */
      const b = Math.min(NB - 1, Math.max(0, Math.floor(r.engagement * NB + 1e-9)));
      counts[r.arm][b]++;
    });
    /* one shared count scale across facets: n is 80 in every arm, so a per-facet scale would
       make the control spike and the designed spread look like the same quantity */
    let top = 1;
    arms.forEach(a => { top = Math.max(top, Math.max.apply(null, counts[a])); });

    const fh = 62, gap = 14, m = {t:34, r:104, b:56, l:104};
    const W = Math.max(560, Math.min(node.clientWidth || 880, 960));
    const H = m.t + arms.length * (fh + gap) - gap + m.b;
    const svg = el('svg', {width:'100%', viewBox:'0 0 ' + W + ' ' + H, role:'img',
      'aria-label':'Distribution of engagement for the designed arm and the three ' +
                   'composition-matched control arms'}, node);
    const x = lin(0, 1, m.l, W - m.r);
    const plotBottom = m.t + arms.length * (fh + gap) - gap;

    /* the filter that does the discriminating, drawn across every facet */
    el('line', {x1:x(ENG_MIN), x2:x(ENG_MIN), y1:m.t - 12, y2:plotBottom + 6, stroke:P.ink2,
                'stroke-width':1, 'stroke-dasharray':'4 3'}, svg);
    txt(svg, x(ENG_MIN) - 7, m.t - 16,
        'engagement filter ≥ ' + ENG_MIN.toFixed(2) + ' — what excludes the trivial signal',
        {'text-anchor':'end', 'font-size':10.5, fill:P.ink2});

    arms.forEach((a, i) => {
      const s = NCSUM[i];
      const yTop = m.t + i * (fh + gap), base = yTop + fh;
      const y = lin(0, top, base, yTop);
      const designed = a === 'designed';
      const col = designed ? P.s1 : P.muted;
      el('line', {x1:m.l, x2:W - m.r, y1:base, y2:base, stroke:P.axis,
                  'shape-rendering':'crispEdges'}, svg);
      grid(svg, m.l, W - m.r, yTop, false);
      txt(svg, m.l - 10, yTop + 9, String(top), {'text-anchor':'end', 'font-size':9.5,
        fill:P.muted});
      txt(svg, m.l - 10, base - 1, '0', {'text-anchor':'end', 'font-size':9.5, fill:P.muted});
      const lab = txt(svg, m.l - 34, yTop + 20, a, {'text-anchor':'end', 'font-size':12,
        fill:designed ? P.ink : P.ink2});
      if (designed) lab.setAttribute('font-weight', '600');
      txt(svg, m.l - 34, yTop + 34, 'n = ' + s.n, {'text-anchor':'end', 'font-size':10,
        fill:P.muted});
      txt(svg, W - m.r + 12, yTop + 20,
          'pass ' + (s.passRate * 100).toFixed(0) + '%', {'font-size':11, fill:P.ink2});
      txt(svg, W - m.r + 12, yTop + 34,
          'median ' + s.engagement.toFixed(2), {'font-size':10, fill:P.muted});

      for (let b = 0; b < NB; b++) {
        const n = counts[a][b];
        const gx = x(b * bw) + 1, gw = Math.max(2, x(bw) - x(0) - 2);
        if (n > 0) {
          const h = base - y(n);
          el('rect', {x:gx, y:y(n), width:gw, height:h, rx:Math.min(3, gw / 2), fill:col}, svg);
        }
        const hit = el('rect', {x:gx - 1, y:yTop, width:gw + 2, height:fh,
                                fill:'transparent'}, svg);
        hit.addEventListener('pointermove', ev => showTip(ev, [
          {value:a, color:col},
          {value:(b * bw).toFixed(2) + ' – ' + ((b + 1) * bw).toFixed(2), name:'engagement'},
          {value:n + ' of ' + s.n + ' constructs'},
          {value:s.note}
        ]));
        hit.addEventListener('pointerleave', hideTip);
      }
      /* median tick, so the summary table's number is visible on the figure too */
      el('line', {x1:x(s.engagement), x2:x(s.engagement), y1:yTop + 4, y2:base,
                  stroke:P.ink, 'stroke-width':2}, svg);
    });
    for (let v = 0; v <= 1.0001; v += 0.1)
      txt(svg, x(v), plotBottom + 18, v.toFixed(1),
          {'text-anchor':'middle', 'font-size':10.5, fill:P.muted});
    txt(svg, (m.l + W - m.r) / 2, H - 18,
        'engagement — expected fraction of DS bases paired to the aptamer',
        {'text-anchor':'middle', 'font-size':11.5, fill:P.ink2});
    txt(svg, m.l - 10, m.t - 16, 'constructs', {'text-anchor':'end', 'font-size':10.5,
      fill:P.ink2});

    const lg = document.createElement('div'); lg.className = 'legend';
    [[P.s1, 'designed displacement strands'],
     [P.muted, 'composition-matched controls (scrambled · reversed · foreign)']]
      .forEach(([c, label]) => {
        const k = document.createElement('span'); k.className = 'k';
        const s = document.createElement('span'); s.className = 'sw'; s.style.background = c;
        k.appendChild(s); k.appendChild(document.createTextNode(label)); lg.appendChild(k);
      });
    const mk = document.createElement('span'); mk.className = 'k';
    const ms = document.createElement('span'); ms.className = 'sw';
    ms.style.width = '3px'; ms.style.background = P.ink;
    mk.appendChild(ms); mk.appendChild(document.createTextNode('arm median'));
    lg.appendChild(mk);
    node.appendChild(lg);
  });
}

function renderControlTable() {
  guard('nc-table', node => {
    if (!NCSUM.length) return empty(node, NC_EMPTY);
    node.textContent = '';
    const head = ['Arm', 'n', 'Median ΔG_switch (kcal/mol)', 'Median closed (fraction)',
                  'Median K_D,app (nM)', 'Median engagement (fraction)',
                  'Pass rate (all filters)', 'Reaching the designed median engagement'];
    const t = document.createElement('table');
    const thead = document.createElement('thead'), tr = document.createElement('tr');
    head.forEach((h, i) => { const th = document.createElement('th');
      th.textContent = h; if (i === 0) th.className = 'l'; tr.appendChild(th); });
    thead.appendChild(tr); t.appendChild(thead);
    const tb = document.createElement('tbody');
    NCSUM.forEach(s => {
      const row = document.createElement('tr');
      const nameTd = document.createElement('td'); nameTd.className = 'l';
      const dot = document.createElement('span');
      dot.className = 'dot'; dot.style.background = s.arm === 'designed' ? P.s1 : P.muted;
      nameTd.appendChild(dot);
      nameTd.appendChild(document.createTextNode(s.arm));
      if (s.arm === 'designed') {
        const b = document.createElement('span'); b.className = 'badge warn';
        b.style.marginLeft = '8px';
        const ic = document.createElement('span'); ic.className = 'ic'; ic.textContent = '⚠';
        b.appendChild(ic);
        b.appendChild(document.createTextNode('pre-selected — see caveat 1'));
        nameTd.appendChild(b);
      }
      nameTd.title = s.note;
      row.appendChild(nameTd);
      [[String(s.n)],
       [s.dg_switch.toFixed(2), 'full range in this arm: ' + s.dgLo.toFixed(2) + ' to ' +
        s.dgHi.toFixed(2) + ' kcal/mol'],
       [s.closed_frac.toFixed(3)],
       [fmt(s.kd_app_nM, 1)],
       [s.engagement.toFixed(2), 'highest in this arm: ' + s.engMax.toFixed(2)],
       [(s.passRate * 100).toFixed(0) + '%'],
       [s.reachRef === null ? '—' : (s.reachRef * 100).toFixed(0) + '%']
      ].forEach(([v, tip]) => { const td = document.createElement('td'); td.textContent = v;
                                if (tip) td.title = tip; row.appendChild(td); });
      tb.appendChild(row);
    });
    t.appendChild(tb); node.appendChild(t);
    const cap = document.createElement('p');
    cap.className = 'sub'; cap.style.marginTop = '8px';
    const d = NCSUM.find(s => s.arm === 'designed');
    const ctl = NCSUM.filter(s => s.arm !== 'designed');
    if (d && ctl.length) {
      const kdLo = Math.min.apply(null, ctl.map(s => s.kd_app_nM));
      const kdHi = Math.max.apply(null, ctl.map(s => s.kd_app_nM));
      cap.textContent =
        'Read K_D,app here with care: the controls look BETTER on it (median ' +
        fmt(kdLo, 0) + '–' + fmt(kdHi, 0) + ' nM against ' + fmt(d.kd_app_nM, 0) +
        ' nM designed) precisely because they barely close, and a switch that never closes ' +
        'never pays the switching penalty. K_D,app is not a discriminator on its own; ' +
        'engagement and closed fraction are, and the pass rate is the joint test. ' +
        'Control engagement never exceeds ' +
        Math.max.apply(null, ctl.map(s => s.engMax)).toFixed(2) +
        ' and no control reaches the designed median of ' + d.engagement.toFixed(2) + '.';
      node.appendChild(cap);
    }
  });
}

/* =========================================================================
   VIEW 11 - cross-model agreement: where each model puts the DNA on IL-6

   Form: a contact track. The question is "do four models cover the same positions on one
   183-residue axis", which is positional-set overlap -- a bar chart of contact counts would
   answer a different question (how many contacts) and hide the entire finding. One lane per
   model, marks where that model places a contact, plus an agreement strip counting how many
   models call each residue. Four models = categorical slots 1-4 in the documented fixed
   order, every lane direct-labelled, so identity never rests on colour alone. Lane order is
   deliberate: the three broadly-converging models sit together and opendde sits last.
   The 4x4 Jaccard matrices are the companion, on one sequential ramp with printed values.
   ========================================================================= */
const CF = DATA.cofold;
const CF_EMPTY = 'aptamer/cofold/cofold_summary.json not yet generated — run ' +
                 'python aptamer/cofold/run_cofold.py, then summarize.py';
const CF_RAMP_MAX = 0.6;   /* fixed Jaccard colour domain so chains stay comparable */
let cfChain = 0;

function cfChainNow() { return (CF && CF.chains[cfChain]) || (CF && CF.chains[0]) || null; }
function cfColor(model) {
  const i = (CF && CF.models ? CF.models : []).indexOf(model);
  return SERIES(P)[(i < 0 ? 0 : i) % 4];
}
function cfPair(map, a, b) {
  const v = map[a + '|' + b];
  return v === undefined ? map[b + '|' + a] : v;
}
const cfNum = (v, d) => (v === null || v === undefined) ? '—' : v.toFixed(d === undefined ? 3 : d);

function renderCofoldTrack() {
  guard('cf-track', node => {
    if (!CF) return empty(node, CF_EMPTY);
    const ch = cfChainNow();
    if (!ch) return empty(node, CF_EMPTY);
    const n = CF.targetLen || 183, off = CF.uniprotOffset;
    const trio = new Set(ch.trio || []), all4 = new Set(ch.consensus || []);
    const counts = ch.counts || {};
    node.textContent = '';

    const cellW = 6, laneH = 19, laneGap = 5, stripH = 34;
    const m = {t:22, r:22, b:104, l:124};
    const lanes = ch.models.length;
    const stripY = m.t + lanes * (laneH + laneGap) + 10;
    const axisY = stripY + stripH + 4;
    const W = m.l + n * cellW + m.r, H = axisY + m.b;
    const svg = el('svg', {width:W, height:H, viewBox:'0 0 ' + W + ' ' + H, role:'img',
      'aria-label':'Contact residues on IL-6 placed by each model for DNA chain ' + ch.dna},
      svg0(node));
    /* 183 residues need the room: scroll the track rather than shrink the cells to nothing */
    svg.style.minWidth = W + 'px';

    /* painted in three passes so nothing hides anything: empty lane beds, then the wash
       columns that mark a three-model call, then the contact marks on top */
    ch.models.forEach((e, i) => {
      el('rect', {x:m.l, y:m.t + i * (laneH + laneGap), width:n * cellW, height:laneH,
                  fill:P.grid, opacity:0.4, rx:4}, svg);
    });
    el('rect', {x:m.l, y:stripY, width:n * cellW, height:stripH, fill:P.grid, opacity:0.4,
                rx:4}, svg);
    (ch.trio || []).forEach(r => {
      el('rect', {x:m.l + (r - 1) * cellW - 0.6, y:m.t - 4, width:cellW + 1.2,
                  height:(axisY - m.t) + 4, fill:P.ink, opacity:0.12, rx:2}, svg);
    });

    ch.models.forEach((e, i) => {
      const y = m.t + i * (laneH + laneGap), c = cfColor(e.model);
      el('rect', {x:4, y:y + laneH / 2 - 5, width:10, height:10, rx:3, fill:c}, svg);
      txt(svg, 20, y + laneH / 2 + 4, e.model, {'font-size':11.5, fill:P.ink});
      txt(svg, m.l - 8, y + laneH / 2 + 4,
          e.n + (e.iptm === null || e.iptm === undefined ? '' : ' · ' + e.iptm.toFixed(2)),
          {'text-anchor':'end', 'font-size':10.5, fill:P.muted});
      e.residues.forEach(r => {
        el('rect', {x:m.l + (r - 1) * cellW + 0.6, y:y + 2, width:cellW - 1.2,
                    height:laneH - 4, rx:2, fill:c}, svg);
      });
    });
    txt(svg, 4, m.t - 8, 'model', {'font-size':10.5, fill:P.muted});
    txt(svg, m.l - 8, m.t - 8, 'contacts · ipTM', {'text-anchor':'end', 'font-size':10.5,
                                                   fill:P.muted});

    /* agreement strip: how many of the four call this residue (sequential, one hue) */
    const maxM = ch.models.length;
    for (let r = 1; r <= n; r++) {
      const k = counts[String(r)] || 0;
      if (!k) continue;
      const h = Math.max(4, (stripH - 6) * k / maxM);
      el('rect', {x:m.l + (r - 1) * cellW + 0.6, y:stripY + stripH - 3 - h, width:cellW - 1.2,
                  height:h, rx:2, fill:seqColor(0.35 + 0.65 * k / maxM)}, svg);
    }
    txt(svg, m.l - 8, stripY + stripH / 2 + 4, 'models agreeing',
        {'text-anchor':'end', 'font-size':10.5, fill:P.ink2});
    txt(svg, m.l - 8, stripY + stripH / 2 + 17, '0–' + maxM,
        {'text-anchor':'end', 'font-size':10, fill:P.muted});

    /* axis: mature-chain numbering, with the UniProt row underneath it. The two numbering
       schemes differ by the signal peptide and misreading them is a real hazard, so both are
       on the axis rather than in a caption. */
    el('line', {x1:m.l, x2:m.l + n * cellW, y1:axisY, y2:axisY, stroke:P.axis,
                'stroke-width':1, 'shape-rendering':'crispEdges'}, svg);
    for (let r = 10; r <= n; r += 10) {
      const x = m.l + (r - 0.5) * cellW;
      el('line', {x1:x, x2:x, y1:axisY, y2:axisY + 4, stroke:P.axis, 'stroke-width':1}, svg);
      txt(svg, x, axisY + 15, String(r), {'text-anchor':'middle', 'font-size':10, fill:P.ink2});
      if (r % 20 === 0)
        txt(svg, x, axisY + 28, String(r + off), {'text-anchor':'middle', 'font-size':9.5,
                                                  fill:P.muted});
    }
    txt(svg, m.l - 8, axisY + 15, 'mature', {'text-anchor':'end', 'font-size':10, fill:P.ink2});
    txt(svg, m.l - 8, axisY + 28, 'UniProt', {'text-anchor':'end', 'font-size':9.5,
                                              fill:P.muted});

    /* selective direct labels: the residues all three converging models call */
    (ch.trio || []).forEach(r => {
      const x = m.l + (r - 0.5) * cellW;
      txt(svg, x, axisY + 42, '▲', {'text-anchor':'middle', 'font-size':8, fill:P.ink});
      const t = txt(svg, 0, 0, String(r), {'font-size':10, fill:P.ink, 'text-anchor':'end'});
      t.setAttribute('transform', 'translate(' + (x + 3.5) + ',' + (axisY + 48) + ') rotate(-90)');
    });

    txt(svg, m.l + n * cellW / 2, H - 26,
        'IL-6 mature chain, residue 1–' + n + ' (1-based within the mature chain)',
        {'text-anchor':'middle', 'font-size':11.5, fill:P.ink2});
    txt(svg, m.l + n * cellW / 2, H - 11,
        'UniProt P05231 numbering = mature residue + ' + off +
        ' (grey row) — the signal peptide 1–' + off + ' is not in these structures',
        {'text-anchor':'middle', 'font-size':10.5, fill:P.muted});

    /* One continuous overlay rather than 183 six-pixel hit columns: at this density a
       per-residue target would be a pinpoint. The pointer picks the nearest residue and a
       crosshair says which one, so there is no dead zone between marks. */
    const cross = el('rect', {x:m.l, y:m.t - 4, width:cellW, height:(axisY - m.t) + 4,
                              fill:P.ink, opacity:0, rx:2}, svg);
    const over = el('rect', {x:m.l, y:m.t - 4, width:n * cellW, height:(axisY - m.t) + 4,
                             fill:'transparent'}, svg);
    over.addEventListener('pointermove', ev => {
      const box = svg.getBoundingClientRect();
      const sx = (ev.clientX - box.left) * (W / (box.width || W));
      const r = Math.max(1, Math.min(n, Math.floor((sx - m.l) / cellW) + 1));
      cross.setAttribute('x', m.l + (r - 1) * cellW);
      cross.setAttribute('opacity', 0.12);
      let aa = '';
      ch.models.forEach(e => { if (e.aa[String(r)]) aa = e.aa[String(r)]; });
      const rows = [{value:'residue ' + r + (aa ? ' ' + aa : ''),
                     name:'UniProt P05231 ' + (r + off)}];
      ch.models.forEach(e => {
        if (e.residues.indexOf(r) >= 0)
          rows.push({color:cfColor(e.model), value:e.model, name:'contact < 4.0 Å'});
      });
      const k = counts[String(r)] || 0;
      rows.push({value:k + ' of ' + maxM + ' models',
                 name:all4.has(r) ? '— every model agrees'
                      : trio.has(r) ? '— ' + CF.trioModels.join(' + ') + ' agree' : ''});
      showTip(ev, rows);
    });
    over.addEventListener('pointerleave', () => {
      cross.setAttribute('opacity', 0); hideTip();
    });

    const lg = document.createElement('div'); lg.className = 'legend';
    const key = (label, style) => {
      const k = document.createElement('span'); k.className = 'k';
      const s = document.createElement('span'); s.className = 'sw';
      Object.assign(s.style, style); k.appendChild(s);
      k.appendChild(document.createTextNode(label)); lg.appendChild(k);
    };
    ch.models.forEach(e => key(e.model + ' contacts', {background:cfColor(e.model)}));
    key('residue all three of ' + CF.trioModels.join(' + ') + ' call (▲ on the axis)',
        {background:'rgba(128,128,128,0.28)'});
    key('agreement strip: taller and darker = more models', {background:seqColor(1)});
    node.appendChild(lg);

    const cap = document.createElement('p');
    cap.className = 'sub'; cap.style.marginTop = '8px';
    cap.textContent =
      ch.dna + ': mean pairwise Jaccard ' + cfNum(ch.meanJ) + ' (' + cfNum(ch.meanJTol) +
      ' allowing ±2 residues) · ' + ch.nConsensus + ' of ' + ch.nUnion +
      ' union residues are called by all ' + maxM + ' models' +
      ((ch.trio || []).length
        ? ' · ' + CF.trioModels.join(' + ') + ' share ' + ch.trio.length + ': ' +
          ch.trio.join(', ')
        : ' · the three converging models share none on this chain') +
      '. ' + (CF.contactNote || '');
    node.appendChild(cap);
  });
}

/* the pairwise matrices: sequential ramp, values printed in every cell (which is also the
   contrast relief), exact on the left and ±2-tolerant on the right as small multiples */
function renderCofoldMatrix() {
  guard('cf-matrix', node => {
    if (!CF) return empty(node, CF_EMPTY);
    const ch = cfChainNow();
    if (!ch) return empty(node, CF_EMPTY);
    node.textContent = '';
    const ms = ch.models.map(e => e.model);
    const cell = 58, lab = 74, gapX = 46;
    const panelW = lab + ms.length * cell;
    const m = {t:44, l:2, b:26};
    const W = m.l + panelW * 2 + gapX, H = m.t + ms.length * cell + m.b;
    const svg = el('svg', {width:W, height:H, viewBox:'0 0 ' + W + ' ' + H, role:'img',
      'aria-label':'Pairwise Jaccard overlap of contact residues between models'}, svg0(node));
    svg.style.minWidth = W + 'px';
    [['exact residue overlap', ch.jaccard, 0], ['tolerant ±2 residues', ch.jaccardTol, 1]]
      .forEach(([title, map, p]) => {
        const x0 = m.l + p * (panelW + gapX);
        txt(svg, x0, 14, title, {'font-size':12, fill:P.ink, 'font-weight':600});
        ms.forEach((b, j) => txt(svg, x0 + lab + j * cell + cell / 2, m.t - 8, b,
                                 {'text-anchor':'middle', 'font-size':10.5, fill:P.ink2}));
        ms.forEach((a, i) => {
          txt(svg, x0 + lab - 8, m.t + i * cell + cell / 2 + 4, a,
              {'text-anchor':'end', 'font-size':10.5, fill:P.ink2});
          ms.forEach((b, j) => {
            const x = x0 + lab + j * cell, y = m.t + i * cell;
            if (i === j) {
              el('rect', {x:x + 1, y:y + 1, width:cell - 2, height:cell - 2, rx:5,
                          fill:P.grid, opacity:0.5}, svg);
              txt(svg, x + cell / 2, y + cell / 2 + 4, '—',
                  {'text-anchor':'middle', 'font-size':12, fill:P.muted});
              return;
            }
            const v = cfPair(map, a, b);
            const fill = (v === undefined || v === null) ? P.grid
                       : seqColor(Math.min(1, v / CF_RAMP_MAX));
            el('rect', {x:x + 1, y:y + 1, width:cell - 2, height:cell - 2, rx:5, fill:fill}, svg);
            const t = txt(svg, x + cell / 2, y + cell / 2 + 4,
                          v === undefined || v === null ? '—' : v.toFixed(3),
                          {'text-anchor':'middle', 'font-size':11.5,
                           fill:(v === undefined || v === null) ? P.muted : inkOn(fill)});
            if (v === 0) t.setAttribute('font-weight', '700');
            const hit = el('rect', {x:x, y:y, width:cell, height:cell, fill:'transparent'}, svg);
            hit.addEventListener('pointermove', ev => showTip(ev, [
              {value:a + ' vs ' + b, name:title},
              {value:v === undefined || v === null ? '—' : v.toFixed(3),
               name:'Jaccard overlap of contact residues' +
                    (v === 0 ? ' — completely disjoint binding sites' : '')}
            ]));
            hit.addEventListener('pointerleave', hideTip);
          });
        });
      });
    const bar = document.createElement('div');
    bar.className = 'rampbar'; bar.style.marginTop = '8px';
    const b = document.createElement('span'); b.className = 'bar';
    b.style.background = 'linear-gradient(90deg,' + seqColor(0) + ',' + seqColor(1) + ')';
    bar.appendChild(document.createTextNode('Jaccard 0'));
    bar.appendChild(b);
    bar.appendChild(document.createTextNode(CF_RAMP_MAX.toFixed(1) +
      ' — same ramp on every chain. 1.000 would be identical contact sets; 0.000 is no shared residue at all.'));
    node.appendChild(bar);
  });
}

function renderCofoldTable() {
  guard('cf-table', node => {
    if (!CF) return empty(node, CF_EMPTY);
    const ch = cfChainNow();
    if (!ch) return empty(node, CF_EMPTY);
    node.textContent = '';
    const t = document.createElement('table');
    const head = ['Model', 'ipTM (interface)', 'pTM (complex)', 'DNA-chain pTM',
                  'Complex pLDDT', 'Contact residues (< 4.0 Å)',
                  'Mean tolerant Jaccard vs all others (all chains)',
                  'Contact patches, mature-chain numbering'];
    const thead = document.createElement('thead'), htr = document.createElement('tr');
    head.forEach((h, i) => { const th = document.createElement('th');
      th.textContent = h; if (i === 0 || i === head.length - 1) th.className = 'l';
      htr.appendChild(th); });
    thead.appendChild(htr); t.appendChild(thead);
    const tb = document.createElement('tbody');
    ch.models.forEach(e => {
      const tr = document.createElement('tr');
      const td0 = document.createElement('td'); td0.className = 'l';
      const dot = document.createElement('span');
      dot.className = 'dot'; dot.style.background = cfColor(e.model);
      td0.appendChild(dot); td0.appendChild(document.createTextNode(e.model));
      tr.appendChild(td0);
      [cfNum(e.iptm), cfNum(e.ptm), cfNum(e.dnaPtm), cfNum(e.plddt), String(e.n),
       cfNum((CF.vsOthers || {})[e.model])].forEach(v => {
        const td = document.createElement('td'); td.textContent = v; tr.appendChild(td); });
      const tdp = document.createElement('td');
      tdp.className = 'l mono'; tdp.textContent = (e.patches || []).join(', ') || '—';
      tdp.style.whiteSpace = 'normal';
      tr.appendChild(tdp);
      tb.appendChild(tr);
    });
    t.appendChild(tb); node.appendChild(t);
    const cap = document.createElement('p');
    cap.className = 'sub'; cap.style.marginTop = '8px';
    const rec = (CF.trioRecurrence || []).filter(r => r[1] > 1);
    cap.textContent =
      'ipTM, pTM and pLDDT are each model’s own confidence and are not comparable across ' +
      'models. The last-but-one column is the only cross-model number here: how much each ' +
      'model’s contact set overlaps the other three, averaged over all ' +
      CF.chains.length + ' DNA chains.' +
      (rec.length ? ' Residues that ' + CF.trioModels.join(' + ') + ' all call on more than ' +
        'one chain: ' + rec.map(r => r[0] + ' (' + r[1] + '/' + CF.chains.length + ' chains)')
        .join(', ') + '. That is a hypothesis to test, not a result.' : '');
    node.appendChild(cap);
  });
}

/* =========================================================================
   VIEW 12 - off-target panels: the metric failing its own control

   Form: emphasis, not categorical. Five targets per panel with one of them the right answer,
   so the on-target bar carries the accent hue, the off-targets are grey, and the dashed rule
   at the on-target value is the thing to read against -- any bar past it is the metric
   ranking a known non-target above the known target. A bar that crosses it takes the reserved
   critical status colour with an icon and a label, never colour alone. Two DNA chains =
   small multiples on one shared x scale.
   ========================================================================= */
function renderOffTarget() {
  guard('ot-chart', node => {
    if (!CF || !(CF.offTarget || []).length)
      return empty(node, CF_EMPTY);
    node.textContent = '';
    const panels = CF.offTarget;
    let maxv = 0;
    panels.forEach(p => p.rows.forEach(r => { maxv = Math.max(maxv, r.iptm || 0); }));
    maxv = Math.max(0.1, maxv) * 1.32;   /* headroom for the value labels */

    const rowH = 32, m = {t:36, r:16, b:44, l:88};
    const panelW = 470;
    const nrows = Math.max.apply(null, panels.map(p => p.rows.length));
    const W = m.l + panelW * panels.length + (panels.length - 1) * 44 + m.r;
    const H = m.t + nrows * rowH + m.b;
    const svg = el('svg', {width:W, height:H, viewBox:'0 0 ' + W + ' ' + H, role:'img',
      'aria-label':'Interface confidence against IL-6 and the gp130-family off-targets'},
      svg0(node));
    svg.style.minWidth = W + 'px';

    panels.forEach((p, pi) => {
      const x0 = m.l + pi * (panelW + 44);
      const x = lin(0, maxv, x0, x0 + panelW - 130);
      const head = p.dna + (p.dnaLen ? ' · ' + p.dnaLen + ' nt' : '') +
                   ' · ' + p.model + ' ipTM';
      txt(svg, x0, 14, head, {'font-size':12.5, fill:P.ink, 'font-weight':600});
      const pkd = (DATA.constructParent || {}).kdNM;
      txt(svg, x0, 28, p.dna.indexOf('parent') === 0
            ? ('the bare parent — the only chain here with a measured K_D' +
               (pkd ? ' (' + pkd.toFixed(0) + ' nM for IL-6)' : ''))
            : 'top-ranked switch construct at the time of the run',
          {'font-size':11, fill:P.muted});

      /* the reference the whole panel is read against */
      const rx = x(p.onTarget || 0);
      el('line', {x1:rx, x2:rx, y1:m.t - 6, y2:m.t + p.rows.length * rowH - 2,
                  stroke:P.ink2, 'stroke-width':1.5, 'stroke-dasharray':'4 3'}, svg);

      p.rows.forEach((r, i) => {
        const y = m.t + i * rowH;
        const beats = !r.onTarget && (r.iptm || 0) > (p.onTarget || 0);
        const fill = r.onTarget ? P.s1 : (beats ? P.crit : P.muted);
        const w = Math.max(2, x(r.iptm || 0) - x0);
        el('rect', {x:x0, y:y + 5, width:w, height:rowH - 14, rx:4, fill:fill}, svg);
        txt(svg, x0 - 8, y + rowH / 2 + 4, r.label,
            {'text-anchor':'end', 'font-size':11.5,
             fill:r.onTarget ? P.ink : (beats ? P.ink : P.ink2),
             'font-weight':(r.onTarget || beats) ? 600 : 400});
        txt(svg, x0 - 8, y + rowH / 2 + 15, r.uniprot || '',
            {'text-anchor':'end', 'font-size':9.5, fill:P.muted});
        let lx = x0 + w + 8;
        const lt = txt(svg, lx, y + rowH / 2 + 4, cfNum(r.iptm),
                       {'font-size':11.5, fill:P.ink});
        lx += 34;
        if (r.onTarget) {
          txt(svg, lx, y + rowH / 2 + 4, 'on-target', {'font-size':10.5, fill:P.ink2});
        } else if (beats) {
          txt(svg, lx, y + rowH / 2 + 4, '⚠ outscores IL-6 by ' +
              Math.abs(r.margin || 0).toFixed(3), {'font-size':10.5, fill:P.crit,
                                                   'font-weight':600});
        } else {
          txt(svg, lx, y + rowH / 2 + 4, (r.margin > 0 ? '+' : '') + cfNum(r.margin),
              {'font-size':10.5, fill:P.muted});
        }
        const hit = el('rect', {x:x0 - 84, y:y, width:panelW, height:rowH,
                                fill:'transparent'}, svg);
        hit.addEventListener('pointermove', ev => showTip(ev, [
          {color:fill, value:r.label + ' (' + (r.uniprot || '?') + ')', name:p.dna},
          {value:'ipTM ' + cfNum(r.iptm),
           name:r.onTarget ? 'the true target' :
                (beats ? 'ranked ABOVE the true target' :
                 'margin ' + (r.margin > 0 ? '+' : '') + cfNum(r.margin) + ' below IL-6')},
          {value:(r.nContacts === null || r.nContacts === undefined) ? '—' : r.nContacts,
           name:'contact residues < 4.0 Å'}
        ]));
        hit.addEventListener('pointerleave', hideTip);
      });
      /* x axis */
      const ay = m.t + p.rows.length * rowH + 2;
      el('line', {x1:x0, x2:x(maxv / 1.32), y1:ay, y2:ay, stroke:P.axis,
                  'stroke-width':1, 'shape-rendering':'crispEdges'}, svg);
      for (let v = 0; v <= maxv / 1.32 + 1e-9; v += 0.2) {
        el('line', {x1:x(v), x2:x(v), y1:ay, y2:ay + 4, stroke:P.axis}, svg);
        txt(svg, x(v), ay + 15, v.toFixed(1), {'text-anchor':'middle', 'font-size':10,
                                               fill:P.muted});
      }
      txt(svg, x0, ay + 30, 'ipTM (0–1, higher = more confident interface) · dashed rule = ' +
          'IL-6, the true target', {'font-size':10.5, fill:P.ink2});
    });

    const lg = document.createElement('div'); lg.className = 'legend';
    [[P.s1, 'IL-6 — the true target'], [P.muted, 'gp130-family off-target, scored below IL-6'],
     [P.crit, '⚠ off-target scored ABOVE IL-6 — the metric failing its control']]
      .forEach(([c, label]) => {
        const k = document.createElement('span'); k.className = 'k';
        const s = document.createElement('span'); s.className = 'sw'; s.style.background = c;
        k.appendChild(s); k.appendChild(document.createTextNode(label)); lg.appendChild(k);
      });
    node.appendChild(lg);
  });
}

function renderOffTargetTable() {
  guard('ot-table', node => {
    if (!CF || !(CF.offTarget || []).length) return empty(node, CF_EMPTY);
    node.textContent = '';
    const panels = CF.offTarget;
    const order = [];
    panels.forEach(p => p.rows.forEach(r => {
      if (order.indexOf(r.target) < 0) order.push(r.target); }));
    const t = document.createElement('table');
    const thead = document.createElement('thead'), htr = document.createElement('tr');
    const head = ['Target', 'UniProt'];
    panels.forEach(p => head.push(p.dna + ' ipTM', p.dna + ' margin vs IL-6',
                                  p.dna + ' contacts (< 4.0 Å)'));
    head.forEach((h, i) => { const th = document.createElement('th');
      th.textContent = h; if (i < 2) th.className = 'l'; htr.appendChild(th); });
    thead.appendChild(htr); t.appendChild(thead);
    const tb = document.createElement('tbody');
    order.forEach(tg => {
      const tr = document.createElement('tr');
      const any = panels.map(p => p.rows.find(r => r.target === tg)).find(Boolean) || {};
      const td0 = document.createElement('td'); td0.className = 'l';
      td0.textContent = any.label || tg;
      if (any.onTarget) {
        const b = document.createElement('span'); b.className = 'badge ok';
        b.style.marginLeft = '8px'; b.textContent = 'on-target';
        td0.appendChild(b);
      }
      tr.appendChild(td0);
      const td1 = document.createElement('td'); td1.className = 'l mono';
      td1.textContent = any.uniprot || '—'; tr.appendChild(td1);
      panels.forEach(p => {
        const r = p.rows.find(q => q.target === tg) || {};
        const beats = r.iptm !== undefined && !r.onTarget && r.iptm > (p.onTarget || 0);
        const a = document.createElement('td'); a.textContent = cfNum(r.iptm);
        if (beats) { a.style.color = P.crit; a.style.fontWeight = '600'; }
        tr.appendChild(a);
        const bcell = document.createElement('td');
        bcell.textContent = r.onTarget ? '—'
          : (r.margin === undefined || r.margin === null ? '—'
             : (r.margin > 0 ? '+' : '') + r.margin.toFixed(3));
        if (beats) { bcell.style.color = P.crit; bcell.style.fontWeight = '600';
                     bcell.title = 'ranked above the true target'; }
        tr.appendChild(bcell);
        const c = document.createElement('td');
        c.textContent = (r.nContacts === undefined || r.nContacts === null) ? '—' : r.nContacts;
        tr.appendChild(c);
      });
      tb.appendChild(tr);
    });
    t.appendChild(tb); node.appendChild(t);
  });
}

function renderCofold() {
  renderCofoldTrack(); renderCofoldMatrix(); renderCofoldTable();
}
function buildCofoldPicker() {
  const sel = document.getElementById('cf-chain');
  if (!sel) return;
  if (!CF) { renderCofold(); renderOffTarget(); renderOffTargetTable(); return; }
  sel.textContent = '';
  CF.chains.forEach((c, i) => {
    const o = document.createElement('option');
    o.value = i;
    o.textContent = c.dna + (c.len ? ' — ' + c.len + ' nt' : '') +
                    (c.kind ? ' ' + c.kind : '') +
                    (c.kdAppNM ? ' · K_D,app ' + c.kdAppNM.toFixed(0) + ' nM at run time' : '');
    sel.appendChild(o);
  });
  sel.value = String(cfChain);
  if (!sel.dataset.wired) {
    sel.addEventListener('change', () => { cfChain = +sel.value || 0; renderCofold(); });
    sel.dataset.wired = '1';
  }
  renderCofold(); renderOffTarget(); renderOffTargetTable();
}

/* ---------- boot + theme ---------- */
function renderAll() {
  P = PAL[mode()];
  renderParents(); renderOccupancy();
  applyFilters();
  /* views 9-10 read their own files and are deliberately NOT scoped by the switch filter row */
  renderMismatchChart(); renderMismatchTable();
  renderControlDist(); renderControlTable();
  renderCofold(); renderOffTarget(); renderOffTargetTable();
  if (window.__ssRedraw) window.__ssRedraw();
  if (window.__molRedraw) window.__molRedraw();
}
function boot() {
  dropFornacGlobalCss();
  buildFilterRow();
  renderParents(); renderOccupancy();
  applyFilters();
  renderMismatchChart(); renderMismatchTable();
  renderControlDist(); renderControlTable();
  buildStructurePicker();
  buildCofoldPicker();
  /* fornac normally executes before this script, but poll briefly in case the CDN is slow;
     if it never arrives the arc diagram that is already on screen simply stays. */
  if (!fornaCtor()) {
    const ft = setInterval(() => {
      if (fornaCtor()) { clearInterval(ft); if (window.__ssRedraw) window.__ssRedraw(); }
    }, 200);
    setTimeout(() => clearInterval(ft), 6000);
  }
  render3D();
  const btn = document.getElementById('theme-toggle');
  btn.addEventListener('click', () => {
    document.documentElement.dataset.theme = mode() === 'dark' ? 'light' : 'dark';
    btn.textContent = mode() === 'dark' ? 'Light mode' : 'Dark mode';
    renderAll();
  });
  btn.textContent = mode() === 'dark' ? 'Light mode' : 'Dark mode';
  matchMedia('(prefers-color-scheme: dark)').addEventListener('change', renderAll);
  let rt; addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(renderAll, 200); });
}
boot();
"""


def stat_tiles(data):
    """Small HTML for the stat row; the hero figure is the passing-construct count."""
    parents = data["parents"] or []
    switches = data["switches"] or []
    tiles = []

    def tile(label, value, unit="", hero=False):
        cls = "value hero" if hero else "value"
        u = f'<span class="unit">{unit}</span>' if unit else ""
        return (f'<div class="stat"><div class="label">{label}</div>'
                f'<div class="{cls}">{value}{u}</div></div>')

    cp = data.get("constructParent")
    mm = data.get("mismatches") or []
    ctl = data.get("controlSummary") or []
    if switches:
        tiles.append(tile("Switch constructs passing the design filters", len(switches), "", hero=True))
        best = min(switches, key=lambda r: r["kd_app_nM"])
        tiles.append(tile("Tightest K_D,app in the library", f"{best['kd_app_nM']:.1f}", "nM"))
        top = switches[0]
        tiles.append(tile("Top-ranked DS", top["ds"], f"window {top['window']}"))
    else:
        tiles.append(tile("Switch constructs", "—", "switches.csv not yet generated"))
    # The scaffold every construct sits on, stated up front: it is a truncation, not one of the
    # full-length rows in view 1.
    if cp:
        unit = (f"{cp['length']} nt — truncation of the {cp['sourceLength']} nt {cp['source']}"
                if cp["truncated"] else f"{cp['length']} nt")
        tiles.append(tile("Parent every construct is built on", cp["name"], unit))
    if mm:
        b = min(mm, key=lambda r: r["kd_app_nM"])
        tiles.append(tile("Tightest single-mismatch variant", f"{b['kd_app_nM']:.1f}",
                          f"nM · {b['name']} at tether {b['tether_nt']} nt"))
    if ctl:
        rates = [f"{s['passRate']:.0%}" for s in ctl if s["arm"] != "designed"]
        tiles.append(tile("Control-arm pass rate", " / ".join(rates) or "—",
                          "scrambled / reversed / foreign"))
    if parents:
        b = min(parents, key=lambda r: r["KD_M"])
        tiles.append(tile("Tightest parent K_D", f"{b['kd_nM']:.1f}", f"nM · {b['name']}"))
        flagged = sum(1 for p in parents if p.get("issues"))
        tiles.append(tile("Parents with uncertain nucleotides", flagged, f"of {len(parents)}"))
    # The two co-folding results, at the top of the page rather than 11 sections down: both are
    # negative, and a reader who only sees the stat row must still see them.
    cf = data.get("cofold")
    if cf:
        worst = min(cf["chains"], key=lambda c: (c["nConsensus"] or 0) / max(1, c["nUnion"] or 1))
        tiles.append(tile("Residues all four models agree the DNA touches",
                          f"{worst['nConsensus']} of {worst['nUnion']}",
                          f"union contacts · {worst['dna']} · view 11"))
        fails = [p for p in cf["offTarget"] if p["beats"]]
        if fails:
            p = fails[0]
            top = p["rows"][0]
            tiles.append(tile("Off-target ranked above the true target",
                              f"{top['label']} {top['iptm']:.3f}",
                              f"vs IL-6 {p['onTarget']:.3f} · {p['dna']} · view 12"))
    return "\n".join(tiles)


def section(title, sub, body_id, extra=""):
    return (f'<section class="card"><header><h2>{title}</h2>'
            f'<p class="sub">{sub}</p></header>{extra}'
            f'<div id="{body_id}"></div></section>')


def build_html(data):
    forna_css = scoped_fornac_css()  # before the payload: it can append a generator note
    # Without the scoped stylesheet fornac would draw strokeless links and its own injected
    # sheet is deleted at runtime, so don't load it at all -- view 7 then uses the arcs.
    forna_script = (
        '<!-- fornac 1.1.8 draws view 7. Note the path: dist/scripts/fornac.js, not\n'
        '     dist/fornac.js (that one 404s). d3 3.5.13 is bundled inside it and jQuery is\n'
        '     never referenced, so it needs no companion script. If this fails to load,\n'
        '     view 7 falls back to the built-in base-pair arc diagram. -->\n'
        '<script src="https://cdn.jsdelivr.net/npm/fornac@1.1.8/dist/scripts/fornac.js">'
        "</script>"
    ) if forna_css else "<!-- fornac not loaded: vendor/fornac.css was unreadable -->"
    payload = json.dumps(data, allow_nan=False).replace("</", "<\\/")
    n_sw = len(data["switches"] or [])
    n_par = len(data["parents"] or [])
    complexes = [s for s in data["structures"] if s["kind"] == "complex"]
    notes = data["notes"]
    note_html = ""
    if notes:
        items = "".join(f"<li>{n}</li>" for n in notes)
        note_html = (f'<div class="issues" style="margin-bottom:20px"><b>Generator notes</b>'
                     f'<ul>{items}</ul></div>')

    complex_note = (
        f"{len(complexes)} co-folded complex(es) found." if complexes else
        "No co-folded aptamer–IL-6 complexes have been generated yet, so only the "
        "IL-6 receptor structure is available here."
    )

    # --- which parent? -------------------------------------------------------------------
    # The constructs sit on a 45-nt truncation while view 1 lists 74-nt full-length sequences,
    # so this banner is stated once at the top of the page and reinforced by the parent table's
    # role column, the switch-library heading and the subtitles of views 9-10.
    cp = data.get("constructParent")
    n_mm = len(data["mismatches"] or [])
    n_ctl = len(data["controls"] or [])
    parent_line = ""
    parent_banner = ""
    if cp and cp["truncated"]:
        kd = f" Its published K_D ({cp['kdNM']:.0f} nM) was measured on this truncated form." \
            if cp.get("kdNM") else ""
        parent_line = (f" · constructs built on <b>{cp['name']}</b> "
                       f"({cp['length']} nt, truncated)")
        parent_banner = (
            f'<div class="callout"><b>Which parent: every construct below is built on '
            f'{cp["name"]}, a {cp["length"]}-nt truncation — not on any sequence in the parent '
            f'table.</b> {cp["name"]} is the first {cp["length"]} nt of the '
            f'{cp["sourceLength"]}-nt <b>{cp["source"]}</b>; module B '
            f'(positions {cp["length"] + 1}–{cp["sourceLength"]}) is removed, which keeps the '
            f'finished construct in the length range real E-AB sensors work at.{kd} '
            f'View 1 shows four full-length parents for reference; the switch library, the '
            f'mismatch variants and the negative controls all use the '
            f'{cp["length"]}-mer only.<br>'
            f'<span class="seq" style="display:inline-block;margin-top:6px">'
            f'{cp["sequence"]}</span> '
            f'<span class="sub">← {cp["name"]}, {cp["length"]} nt, 5′→3′</span></div>')
    elif cp:
        parent_line = f" · constructs built on <b>{cp['name']}</b> ({cp['length']} nt)"
    switch_head = (f' — every construct built on {cp["name"]}, {cp["length"]} nt'
                   + (f' (truncation of {cp["source"]})' if cp and cp["truncated"] else "")
                   ) if cp else ""

    # --- view 9 headline ------------------------------------------------------------------
    mm_call = ""
    mm_parent = (f"{n_mm} variants" if n_mm else "none generated yet")
    if n_mm and cp:
        teth = sorted({r["tether_nt"] for r in data["mismatches"]})
        mm_parent = (f"{n_mm} variants at tether {teth[0]}–{teth[-1]} nt, all on "
                     f"{cp['name']} ({cp['length']} nt"
                     + (f", truncation of {cp['source']})" if cp["truncated"] else ")"))
    h = data.get("mismatchHeadline")
    if h:
        bc, lb, ma = h["bestCovered"], h["libBest"], h["matchAlt"]
        equiv = (f"Matching that {bc['kd_app_nM']:.0f} nM by lengthening instead — with at "
                 f"least the same {bc['rand_covered']} randomised nucleotides covered — means "
                 f"a <b>{ma['tether_nt']}-nt tether</b> (DS <code>{ma['ds']}</code> at window "
                 f"{ma['window']}, {ma['kd_app_nM']:.0f} nM), a "
                 f"{ma['tether_nt'] - bc['tether_nt']}-nt jump"
                 if ma else
                 "No unmutated construct covering as much reaches that K_D,app at any tether")
        lib_bit = (f"; the tightest such construct anywhere in the library "
                   f"({lb['kd_app_nM']:.0f} nM, DS <code>{lb['ds']}</code> at window "
                   f"{lb['window']}) pays a <b>{lb['tether_nt']}-nt tether</b>." if lb else ".")
        mm_call = (
            f'<div class="callout"><b>{h["improved"]} of {h["n"]} passing variants tighten '
            f'K_D,app, and none of them lengthens the tether.</b> '
            f'{h["n"]} single-mismatch variants of {h["nBase"]} short-tether constructs '
            f'(tether {h["tetherLo"]}–{h["tetherHi"]} nt). Best overall: '
            f'<code>{h["best"]["name"]}</code>, '
            f'{h["best"]["kd_app_wt_nM"]:.0f} → <b>{h["best"]["kd_app_nM"]:.0f} nM</b> at '
            f'tether {h["best"]["tether_nt"]} nt. Best of those covering '
            f'{h["covered"]} randomised nt: <code>{bc["name"]}</code>, '
            f'{bc["kd_app_wt_nM"]:.0f} → <b>{bc["kd_app_nM"]:.0f} nM</b> at tether '
            f'{bc["tether_nt"]} nt. {equiv}{lib_bit} '
            f'A longer tether is not free: it slows reclosure and the response, which is why '
            f'the same affinity bought with a mismatch is worth more than the same affinity '
            f'bought with linker.</div>')

    # --- view 10 caveats ------------------------------------------------------------------
    # Static HTML, above the figure, and not a footnote: this is validation evidence and both
    # of these limit what it proves. They stay legible even if the chart JS fails.
    nc_caveats = ""
    summ = {s["arm"]: s for s in (data.get("controlSummary") or [])}
    des = summ.get("designed")
    ctls = [s for s in (data.get("controlSummary") or []) if s["arm"] != "designed"]
    nc_lead = ("Designed displacement strands and their controls"
               if not des else
               f"{des['n']} designed displacement strands, sampled at a fixed stride across "
               f"the whole ranked library so the sample is not just its top"
               + (f" (of {n_sw} passing constructs)" if n_sw else ""))
    if des and ctls:
        rates = ", ".join(f"{s['arm']} {s['passRate']:.0%}" for s in ctls)
        rate_lo = min(s["passRate"] for s in ctls)
        rate_hi = max(s["passRate"] for s in ctls)
        dg_lo = min(s["dg_switch"] for s in ctls)
        dg_hi = max(s["dg_switch"] for s in ctls)
        eng_max = max(s["engMax"] for s in ctls)
        nc_caveats = (
            f'<div class="caveat"><b>⚠ Caveat 1 — the designed arm\'s '
            f'{des["passRate"]:.0%} pass rate is partly circular.</b> Those constructs were '
            f'sampled from switches.csv, which is the set that already passed exactly these '
            f'filters, so of course they pass again; that number measures nothing. The '
            f'non-circular evidence is the rest of the figure: the <b>engagement gap</b> '
            f'(designed median {des["engagement"]:.2f} against '
            f'{min(s["engagement"] for s in ctls):.2f}–'
            f'{max(s["engagement"] for s in ctls):.2f} for the controls, with no control '
            f'exceeding {eng_max:.2f}) and the <b>control pass rates of '
            f'{rate_lo:.0%}–{rate_hi:.0%}</b> ({rates}) — the controls were never '
            f'pre-selected for anything, so their failure is a real result and the designed '
            f'arm\'s success is not.</div>'
            f'<div class="caveat"><b>⚠ Caveat 2 — the controls still show a negative '
            f'ΔG_switch, {dg_lo:+.2f} to {dg_hi:+.2f} kcal/mol at the arm medians.</b> '
            f'Appending any DNA lowers a construct\'s free energy, because the tail forms '
            f'<i>some</i> structure whether or not it is complementary to the aptamer. That '
            f'residual is exactly the trivial signal that raw ΔG_switch on its own would '
            f'admit, and the <b>engagement filter (≥ {BUDGET.get("MIN_ENGAGEMENT", 0.6):.2f}, '
            f'the dashed rule)</b> is what excludes it. Read the '
            f'{rate_lo:.0%}–{rate_hi:.0%} control pass rate as this pipeline\'s '
            f'false-positive floor, not as zero.</div>')

    # --- views 11-12: co-folding ----------------------------------------------------------
    # Every caveat below is static HTML above the figure it qualifies, with its numbers read
    # from cofold_summary.json rather than typed in, so a re-run of the co-folding step moves
    # the caveat too. None of this is a footnote: each one changes how the figure reads.
    cf = data.get("cofold")
    stale = data.get("cofoldStale")
    cofold_head = "Co-folded structures — not yet generated"
    cofold_caveats = cofold_lead = cf_caveat = ot_headline = ot_caveat = ""
    cofold_n = 0
    cf_model_count = "the models"
    if cf:
        cf_model_count = {1: "one model", 2: "two models", 3: "three models",
                          4: "four models"}.get(len(cf["models"]),
                                                f"{len(cf['models'])} models")
        c = cf["counts"]
        cofold_n = c.get("succeeded", 0)
        cofold_head = (
            f"Co-folded structures — {cofold_n} of "
            f"{c.get('total_predictions_attempted', cofold_n)} predictions, "
            f"{len(cf['models'])} models, {len(cf['chains'])} DNA chains")
        cofold_lead = (
            '<div class="callout"><b>These are figures and a hypothesis. They are not ranking '
            'input, and nothing below feeds the pipeline\'s ranking.</b> '
            'The ranking comes from switches.csv — thermodynamics and a measured K_D — and '
            'these predictions were run after it, never into it. Read views 11 and 12 as two '
            'negative results: the models do not agree with each other, and the one '
            'specificity check with a known right answer gets it wrong.'
            + (f'<br><span class="sub" style="display:inline-block;margin-top:6px">'
               f'From <code>cofold_summary.json</code>: {cf["interpretation"]}</span>'
               if cf.get("interpretation") else "")
            + '</div>')

        # 1. the constructs are stale -----------------------------------------------------
        if stale and stale["now"]:
            t, now, sh = stale["then"], stale["now"], stale["short"]
            nowlen = f"{now['len']} nt" if now["len"] else "a shorter construct"
            shorter = (f' and <b>{sh["len"]} nt at {sh["kd"]:.0f} nM</b> ({sh["what"]}) — '
                       f'{t["len"] - sh["len"]} nt shorter than {t["dna"]} and still tighter '
                       f'than every chain co-folded here'
                       if sh and t.get("len") and sh["len"] < t["len"] else "")
            cofold_caveats += (
                f'<div class="caveat"><b>⚠ The co-folded constructs are stale — do not compare '
                f'them to views 3–9.</b> They were snapshotted before the ensemble-ΔG and '
                f'K_closed corrections and before mismatch tuning existed. '
                f'<b>{t["dna"]}</b>, the top-ranked construct then, is '
                f'{t["len"]} nt at K_D,app {t["kd"]:.0f} nM; the {stale["nUsed"]} constructs '
                f'here span {stale["thenLenLo"]}–{stale["thenLenHi"]} nt and '
                f'{stale["thenKdLo"]:.0f}–{stale["thenKdHi"]:.0f} nM. The library now reaches '
                f'<b>{nowlen} at {now["kd"]:.0f} nM</b> ({now["what"]}){shorter}. Same '
                f'pipeline, different molecules: every structure below is of a sequence the '
                f'shortlist has since moved past.</div>')
        # 2. opendde's flat confidence ----------------------------------------------------
        od = [e["iptm"] for ch in cf["chains"] for e in ch["models"]
              if e["model"] == "opendde" and e["iptm"] is not None]
        vs = cf["vsOthers"] or {}
        if od and len(od) > 1:
            cofold_caveats += (
                f'<div class="caveat"><b>⚠ opendde returns a near-constant ipTM '
                f'({min(od):.2f}–{max(od):.2f}) whatever you feed it.</b> Across all '
                f'{len(od)} DNA chains its interface confidence moves by '
                f'{max(od) - min(od):.2f}, while it is also the model that agrees least with '
                f'the others'
                + (f' (mean tolerant Jaccard {vs["opendde"]:.3f} against the other three)'
                   if "opendde" in vs else "")
                + '. It is the most confident model and the least concordant one. A '
                  'confidence score that does not vary with its input is not carrying '
                  'information about the input — read it as a warning sign, not a strength.'
                  '</div>')
        # 3. ipTM collapses with tether length --------------------------------------------
        bz = {ch["dna"]: (e["iptm"], ch["len"]) for ch in cf["chains"] for e in ch["models"]
              if e["model"] == "boltz2" and e["iptm"] is not None}
        par = next((v for k, v in bz.items() if k.startswith("parent")), None)
        con = [v for k, v in bz.items() if not k.startswith("parent")]
        if par and con:
            lo = min(v[0] for v in con)
            hi = max(v[0] for v in con)
            cofold_caveats += (
                f'<div class="caveat"><b>⚠ ipTM collapses with tether length, which is most '
                f'likely an artefact rather than worse binding.</b> Boltz-2 scores '
                f'{par[0]:.3f} on the bare {par[1]}-nt parent and {lo:.3f}–{hi:.3f} on the '
                f'{min(v[1] for v in con)}–{max(v[1] for v in con)}-nt constructs. The '
                f'difference is mostly {min(v[1] for v in con) - par[1]}–'
                f'{max(v[1] for v in con) - par[1]} nt of single-stranded poly-T spacer being '
                f'scored: it has no defined conformation, so whatever the model draws there is '
                f'arbitrary and it drags the interface score down by existing. Compare '
                f'constructs with constructs; never parent against construct.</div>')

        # view 11's own caveat -------------------------------------------------------------
        zeros = sorted({tuple(sorted(k.split("|"))) for ch in cf["chains"]
                        for k, v in ch["jaccard"].items() if v == 0})
        jl = [ch["meanJ"] for ch in cf["chains"] if ch["meanJ"] is not None]
        jt = [ch["meanJTol"] for ch in cf["chains"] if ch["meanJTol"] is not None]
        agree_max = max((ch["nConsensus"] or 0) for ch in cf["chains"])
        if jl:
            pairs = ", ".join(a + " vs " + b for a, b in zeros[:3])
            cf_caveat = (
                f'<div class="caveat"><b>⚠ The four models do not agree on where the aptamer '
                f'binds, so the binding mode is unresolved.</b> Mean pairwise Jaccard of the '
                f'contact residues is {min(jl):.3f}–{max(jl):.3f} '
                f'({min(jt):.3f}–{max(jt):.3f} allowing ±2 residues), and at most '
                f'{agree_max} residue out of unions of '
                f'{min(ch["nUnion"] for ch in cf["chains"])}–'
                f'{max(ch["nUnion"] for ch in cf["chains"])} is called by all four'
                + (f'. Some pairs share nothing at all: {pairs} score exactly 0.000.'
                   if zeros else '.')
                + ' Three of them — ' + " + ".join(cf["trioModels"]) + ' — do broadly converge '
                  'on two patches, and those recurring residues are the obvious thing to probe '
                  'experimentally, but three models that share training data and biases are '
                  'correlated evidence, not independent replication.</div>')

        # view 12's headline: this one must be impossible to miss --------------------------
        fails = [p for p in cf["offTarget"] if p["beats"]]
        tight = None
        for p in cf["offTarget"]:
            for r in p["rows"]:
                if not r["onTarget"] and r["margin"] is not None and r["margin"] > 0:
                    if tight is None or r["margin"] < tight[0]:
                        tight = (r["margin"], r["label"], p["dna"])
        if fails:
            p = fails[0]
            top = max((r for r in p["rows"] if not r["onTarget"]),
                      key=lambda r: r["iptm"] or 0)
            # p["rows"] is already sorted by ipTM, so IL-6's index is its rank among the five
            place = [r["target"] for r in p["rows"]].index("IL6") + 1
            place = {1: "first", 2: "second", 3: "third", 4: "fourth",
                     5: "fifth"}.get(place, f"#{place}")
            ot_headline = (
                f'<div class="caveat"><b>⚠ The metric ranks a known non-target above the known '
                f'target.</b> For <b>{p["dna"]}</b> — the only chain here with a real measured '
                + (f'affinity (K_D {cp["kdNM"]:.0f} nM for IL-6 by SPR, and no measurable '
                   f'binding to human serum albumin) — ' if cp and cp.get("kdNM")
                   else "affinity — ")
                + f'{p["model"]} scores '
                f'<b>{top["label"]} {top["iptm"]:.3f}</b> against <b>IL-6 '
                f'{p["onTarget"]:.3f}</b>. The true target comes {place} of the '
                f'{len(p["rows"])} cytokines, behind one the aptamer is not known to bind, by '
                f'{abs(top["margin"]):.3f} ipTM. This is a negative result <i>about the '
                f'method</i>, not specificity data'
                + (f': because the metric fails the one case where the right answer is known, '
                   f'the {tight[1]} margin of +{tight[0]:.3f} on {tight[2]} carries no '
                   f'information either — it is well inside the error the control just '
                   f'demonstrated.' if tight else '.')
                + '</div>')
        ot_caveat = (
            '<div class="caveat"><b>⚠ Do not quote any number in this view as specificity.</b> '
            'The design is the most defensible thing here — same model, same DNA, same '
            'protocol, five homologous cytokines, so systematic errors partly cancel — and it '
            'still fails. What the panels rule out is the metric as a selection filter; they '
            'are not evidence that the aptamer is or is not selective. The two panels do not '
            'even agree with each other on the ordering of the off-targets.</div>')

    # --- view 8: which complexes are inlined, and what is not ------------------------------
    inlined = [s for s in data["structures"] if s.get("chainA")]
    on_disk = (cf or {}).get("nOnDisk", 0)
    if inlined:
        left = max(0, on_disk - len(inlined))
        complex_note = (
            f"{len(inlined)} co-folded complexes are inlined here, chosen because they carry "
            f"the argument: the same DNA chain placed by all four models (view 11's "
            f"disagreement, in 3D) plus the off-target run that outscores the true target "
            f"(view 12). Chain A is the protein and chain B the DNA, coloured apart. "
            + (f"The other {left} of the {on_disk} predicted structures are not inlined — each "
               f".cif is ~250 KB and all of them would add several MB to this page — they are "
               f"on disk in <code>aptamer/cofold/structures/</code>."
               if left else "")
            + " These are predictions, not solved structures; see the caveats in views 11–12.")
    elif complexes:
        complex_note = f"{len(complexes)} co-folded complex(es) found."

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>IL-6 aptamer switch design — dashboard</title>
<style data-dash>{CSS}
/* vendor/fornac.css, every selector rewritten to sit under {FORNA_SCOPE} by
   dashboard.py:scoped_fornac_css(). Unscoped it would resize every chart on the page. */
{forna_css}</style>
</head>
<body data-palette="#2a78d6,#eb6834,#1baf7a,#eda100">
<div class="viz-root">

<header class="top">
  <div class="grow">
    <h1>Structure-switching IL-6 aptamers</h1>
    <p class="sub">{n_par} parent aptamers · {n_sw} ranked ISD switch constructs ·
       {n_mm} single-mismatch variants · {n_ctl} negative-control rows ·
       {cofold_n} co-folded complexes{parent_line} ·
       generated by <code>aptamer/dashboard.py</code>. All data is inlined — works offline
       except the 3D panel.</p>
  </div>
  <button id="theme-toggle" type="button">Dark mode</button>
</header>

{note_html}

{parent_banner}

<div class="stats">{stat_tiles(data)}</div>

{section("1 · Parent aptamers",
         "Published Neomer candidates, reconstructed and folded — full-length sequences, none of "
         "which is the construct scaffold. The role column says which row the constructs' "
         "truncated parent comes from, and the struck-through nucleotides are the module the "
         "truncation drops. Reconstruction warnings flag a genuinely uncertain nucleotide — "
         "those constructs must not be ordered blind.",
         "parents-body")}

<section class="card">
  <header><h2>2 · Fractional occupancy in wound fluid</h2>
  <p class="sub">Fraction of aptamer bound vs IL-6 concentration, from each published K_D.
  Vertical rules mark the wound-fluid medians; the wash spans post-surgical to infected;
  the amber rule is the skin-model reading (5 ng/mL).</p></header>
  <div id="occ-chart"></div>
  <details class="tv" open><summary>Table view — occupancy at each landmark</summary>
    <div id="occ-table" style="margin-top:8px"></div></details>
</section>

<h2 style="margin:28px 0 10px">Switch library{switch_head}</h2>
<p class="sub" style="margin-bottom:12px">These filters scope views 3–7 — the table, the
scatter, the heatmap, the coverage track and the switch structures all re-render against the
same slice. Views 9 and 10 read their own files and are not filtered here.</p>
<div class="filters" id="sw-filters"></div>

{section("3 · Ranked switch constructs",
         "Default order is the pipeline's ranking: randomised nucleotides covered, then "
         "engagement, then short DS. Click any header to sort.",
         "sw-table")}

<div class="grid2">
  {section("4 · Gain vs affinity",
           "Every construct: stabilisation energy against the affinity it costs. Colour = "
           "randomised nucleotides the DS covers; size = engagement.",
           "sw-scatter")}
  {section("5 · Linker × DS length",
           "Best (lowest) K_D,app reachable at each combination — the two knobs the "
           "design turns.",
           "sw-heatmap")}
</div>

{section("6 · Where the DS windows land",
         "The parent sequence as an axis. Bars count the passing constructs whose displacement "
         "strand targets each position; the amber wash marks the randomised library blocks.",
         "sw-coverage")}

<section class="card">
  <header><h2>7 · Secondary structure</h2>
  <p class="sub">The real 2D layout, drawn by <a href="https://github.com/ViennaRNA/fornac"
  target="_blank" rel="noopener">fornac</a> from the dot-bracket the pipeline folded —
  stems, loops and junctions where they actually sit. Nothing is re-folded here and the DNA
  sequence is shown verbatim (no T→U conversion). Nucleotides are coloured by construct
  region: aptamer, poly-T linker, displacement strand. Drag to pan, scroll to zoom.
  fornac comes from a CDN; without it this panel falls back to the built-in base-pair arc
  diagram and says so, so it still works offline.</p></header>
  <div style="margin-bottom:10px">
    <label class="sub" for="ss-select">Structure&nbsp;</label>
    <select id="ss-select" style="max-width:min(560px,100%)"></select>
  </div>
  <div id="ss-chart"></div>
</section>

<section class="card">
  <header><h2>8 · 3D structures</h2>
  <p class="sub">{complex_note} This panel loads 3Dmol.js from a CDN, so it is the one view
  that needs internet.</p></header>
  <div style="margin-bottom:10px">
    <label class="sub" for="mol-select">Structure&nbsp;</label>
    <select id="mol-select" style="max-width:min(620px,100%)"></select>
  </div>
  <div class="viewer" id="mol-viewer"></div>
  <div class="legend" id="mol-key"></div>
</section>

<section class="card">
  <header><h2>9 · Mismatch refinement — affinity that costs no tether</h2>
  <p class="sub">Single-mismatch variants of the shortest-tether constructs
  ({mm_parent}). Each row is one variant: the open mark is the unmutated construct's K_D,app,
  the filled mark is the same construct with one base changed in the displacement strand, and
  the tether column does not move between them. DS length and linker length both trade
  affinity against response speed; a mismatch weakens the duplex while leaving both lengths
  alone, so it is the only knob that improves K_D,app without slowing the sensor. Values are
  read from <code>mismatches.csv</code> as the pipeline produced them. Not scoped by the
  filter row above.</p></header>
  {mm_call}
  <div id="mm-chart"></div>
  <details class="tv" open><summary>Table view — every passing variant, sortable by any
    column</summary>
    <div id="mm-table" style="margin-top:8px"></div></details>
</section>

<section class="card">
  <header><h2>10 · Negative controls — does the score detect complementarity or just DNA?</h2>
  <p class="sub">{nc_lead}, each against three controls matched to it on
  length and base composition: <b>scrambled</b> (same bases shuffled), <b>reversed</b> (same
  bases, wrong pairing register) and <b>foreign</b> (the reverse complement of a window of a
  <i>shuffled</i> aptamer — a genuine duplex-former that simply is not complementary to this
  aptamer, and so the strict control). One facet per arm, shared axis; the designed arm is in
  the accent hue and the controls in grey because the gap between them is the result.
  Not scoped by the filter row above.</p></header>
  {nc_caveats}
  <div id="nc-chart"></div>
  <details class="tv" open><summary>Table view — median score and pass rate per arm</summary>
    <div id="nc-table" style="margin-top:8px"></div></details>
</section>

<h2 style="margin:28px 0 10px">{cofold_head}</h2>
<p class="sub" style="margin-bottom:12px">Read from <code>aptamer/cofold/</code>, which this
page never writes. Views 11 and 12 are not scoped by the switch filter row: the co-folding step
ran against its own snapshot of the library.</p>

{cofold_lead}
{cofold_caveats}

<section class="card">
  <header><h2>11 · Cross-model agreement — where {cf_model_count} put the DNA on IL-6</h2>
  <p class="sub">One lane per structure-prediction model, along the IL-6 mature chain as a
  residue axis. A mark is an IL-6 residue with at least one heavy atom within 4.0 Å of any DNA
  atom in that model's predicted complex — measured from the returned coordinates, not a score
  the model reports. The strip underneath counts how many models call each residue, and ▲ marks
  the residues all three of the converging models agree on. The matrices are the same contact
  sets as pairwise Jaccard overlap, exact and allowing ±2 residues (two models can hit the same
  patch and still score near zero if their lists are offset by a residue or two).</p></header>
  {cf_caveat}
  <div style="margin-bottom:10px">
    <label class="sub" for="cf-chain">DNA chain&nbsp;</label>
    <select id="cf-chain" style="max-width:min(620px,100%)"></select>
  </div>
  <div id="cf-track"></div>
  <div id="cf-matrix" style="margin-top:18px"></div>
  <details class="tv" open><summary>Table view — per-model confidence, contact count and
    patches</summary>
    <div id="cf-table" style="margin-top:8px"></div></details>
</section>

<section class="card">
  <header><h2>12 · Off-target panels — the specificity metric fails its own control</h2>
  <p class="sub">Boltz-2 ipTM for one DNA chain against IL-6 and the four other gp130-family
  cytokines, same model and same protocol throughout. Two chains, two panels, one shared
  scale. The dashed rule in each panel is the on-target IL-6 value; an off-target bar that
  crosses it is the model preferring a cytokine the aptamer is not known to bind.</p></header>
  {ot_headline}
  {ot_caveat}
  <div id="ot-chart"></div>
  <details class="tv" open><summary>Table view — ipTM, margin and contact count for every
    target</summary>
    <div id="ot-table" style="margin-top:8px"></div></details>
</section>

<footer class="foot">
Regenerate with <code>python aptamer/dashboard.py</code>. Never hand-edit dashboard.html.
Values are read as the pipeline produced them; only residence time (1/k_off) and fractional
occupancy are derived for display.
</footer>
</div>

<script id="dash-data" type="application/json">{payload}</script>
<script src="https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"></script>
{forna_script}
<script>{JS}</script>
</body>
</html>
"""


def main():
    data = build_data()
    html = build_html(data)
    OUT.write_text(html, encoding="utf-8")
    kb = OUT.stat().st_size / 1024
    cp = data.get("constructParent")
    print(f"parents.json   : {len(data['parents'] or [])} aptamers")
    print(f"switches.csv   : {len(data['switches'] or [])} constructs")
    print(f"mismatches.csv : {len(data['mismatches'] or [])} variants"
          if data["mismatches"] is not None else
          "mismatches.csv : not present -> view 9 shows 'not yet generated'")
    print(f"neg controls   : {len(data['controls'] or [])} rows, "
          f"{len(data['controlSummary'])} arms"
          if data["controls"] is not None else
          "neg controls   : not present -> view 10 shows 'not yet generated'")
    if cp:
        print(f"construct parent: {cp['name']} ({cp['length']} nt"
              + (f", truncation of {cp['source']} at {cp['sourceLength']} nt)"
                 if cp["truncated"] else ")"))
    cf = data.get("cofold")
    if cf:
        c = cf["counts"]
        print(f"cofold summary : {c.get('succeeded')}/{c.get('total_predictions_attempted')} "
              f"predictions, {len(cf['models'])} models, {len(cf['chains'])} DNA chains, "
              f"{len(cf['offTarget'])} off-target panel(s)")
        for ch in cf["chains"]:
            print(f"  {ch['dna']:<11}: mean Jaccard {ch['meanJ']}, "
                  f"{ch['nConsensus']}/{ch['nUnion']} residues all models agree on, "
                  f"trio consensus {ch['trio'] or '-'}")
        for p in cf["offTarget"]:
            print(f"  {p['dna']:<11}: on-target {p['onTarget']:.3f}, min margin "
                  f"{p['minMargin']}" + (f", OUTSCORED BY {', '.join(p['beats'])}"
                                         if p["beats"] else ""))
    else:
        print("cofold         : aptamer/cofold/cofold_summary.json absent -> views 11-12 show "
              "'not yet generated'")
    for s in data["structures"]:
        print(f"{s['name']:<28}: {s['atoms']} atoms ({s['kind']})")
    if not any(s["kind"] == "complex" for s in data["structures"]):
        print("co-folded      : none yet -> 3D panel shows the receptor only")
    for n in data["notes"]:
        print(f"note           : {n}")
    print(f"\nwrote {OUT} ({kb:.0f} KB)")


if __name__ == "__main__":
    main()
