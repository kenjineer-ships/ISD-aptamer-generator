"""Fine-tune top switch constructs with single mismatches in the displacement strand.

The third design axis. DS length and linker length both change the *kinetics* while they change
the thermodynamics -- a longer DS is slower to displace, a longer linker slows reclosure. A
mismatch inside the DS weakens the duplex while leaving both lengths alone, so it is the only
knob that moves dG_switch without moving the response speed. Wilson 2019 names these variants
`9_33 p7CG`, i.e. DS 9, loop 33, position 7 C->G; the same convention is used here.

Applied as refinement to already-good candidates rather than as a fourth sweep axis: a global
enumeration would be ~28x the grid for a knob that only matters once a design is chosen.

Reads  aptamer/switches.csv  (already rank-ordered)
Writes aptamer/mismatches.csv
"""

import csv
import math
import pathlib

import RNA
from build_parents import fold, parent
from switch_library import (
    CLOSED_MIN,
    KD_APP_MAX_NM,
    MIN_ENGAGEMENT,
    PARENT,
    RT,
    SELECTIVITY_MIN,
    dimer_screen,
    engagement,
    k_closed_from,
    max_homopolymer,
    revcomp,
)

TOP_N = 40  # how many ranked constructs to refine
HERE = pathlib.Path(__file__).parent


def score(aptamer, g_aptamer, ds, window_start, linker_len):
    """Everything switch_library computes, for an arbitrary (possibly mismatched) DS."""
    n = len(aptamer)
    tail = "T" * linker_len + ds
    _, _, g_closed, fc = fold(aptamer + tail)
    g_open = g_aptamer + fold(tail)[2]
    dg = g_closed - g_open
    k_closed = k_closed_from(dg)  # imported, not reimplemented: the -1 must not drift
    designed, homo, hairpin = dimer_screen(ds, aptamer[window_start:window_start + len(ds)])
    return dict(
        dg_switch=round(dg, 2),
        closed_frac=round(k_closed / (1 + k_closed), 3),
        kd_app_nM=round(KD_PARENT * (1 + k_closed) * 1e9, 1),
        engagement=round(engagement(fc, n, n + linker_len, len(ds)), 2),
        selectivity=round(min(homo, hairpin) - designed, 2),
        tether_nt=linker_len + len(ds),
    )


def passes(s, ds):
    return (s["closed_frac"] >= CLOSED_MIN
            and s["kd_app_nM"] <= KD_APP_MAX_NM
            and s["engagement"] >= MIN_ENGAGEMENT
            and s["selectivity"] >= SELECTIVITY_MIN
            and max_homopolymer(ds) <= 4)


RNA.params_load_DNA_Mathews2004()
APTAMER, KD_PARENT = parent(PARENT)
G_APTAMER = fold(APTAMER)[2]


def main():
    src = HERE / "switches.csv"
    all_rows = list(csv.DictReader(src.open()))
    # Refine the SHORTEST-TETHER constructs, not the top-ranked ones. The rank leaders are
    # already at the closed-fraction floor, so weakening their DS drops them out of the filter
    # entirely -- mismatches cannot help a design that is already as weak as allowed. Short
    # tethers are the opposite case: fast and compact, but strongly closed and therefore poor
    # on affinity. A mismatch buys back affinity without touching either length.
    base_rows = sorted(all_rows, key=lambda r: int(r["tether_nt"]))[:TOP_N]
    print(f"{PARENT}: refining the {len(base_rows)} shortest-tether constructs of {src.name} "
          f"(tether {base_rows[0]['tether_nt']}-{base_rows[-1]['tether_nt']} nt) "
          f"with single DS mismatches\n")

    out, improved = [], 0
    for base in base_rows:
        ds_wt = base["ds"]
        start = int(base["window"].split("-")[0])
        linker_len = int(base["linker_len"])
        wt = score(APTAMER, G_APTAMER, ds_wt, start, linker_len)

        best = None
        for pos in range(1, len(ds_wt) + 1):          # 1-indexed, Wilson's convention
            for alt in "ACGT":
                if alt == ds_wt[pos - 1]:
                    continue
                ds = ds_wt[:pos - 1] + alt + ds_wt[pos:]
                s = score(APTAMER, G_APTAMER, ds, start, linker_len)
                if not passes(s, ds):
                    continue
                rec = dict(name=f"{len(ds_wt)}_{linker_len}_p{pos}{ds_wt[pos - 1]}{alt}",
                           ds_wt=ds_wt, ds=ds, window=base["window"],
                           linker_len=linker_len, mismatch_pos=pos,
                           kd_app_wt_nM=wt["kd_app_nM"],
                           d_kd_nM=round(s["kd_app_nM"] - wt["kd_app_nM"], 1),
                           rand_covered=int(base["rand_covered"]), **s)
                out.append(rec)
                if best is None or rec["kd_app_nM"] < best["kd_app_nM"]:
                    best = rec
        if best and best["kd_app_nM"] < wt["kd_app_nM"]:
            improved += 1

    out.sort(key=lambda r: (-r["rand_covered"], r["kd_app_nM"], r["tether_nt"]))
    print(f"{len(out)} mismatch variants pass all filters; "
          f"{improved}/{len(base_rows)} parent constructs improved by their best mismatch\n")
    print(f"{'variant':<18}{'DS':<13}{'dG':>7}{'closed':>8}{'KDapp':>8}"
          f"{'dKD':>8}{'eng':>6}{'selec':>7}{'teth':>6}")
    for r in out[:15]:
        print(f"{r['name']:<18}{r['ds']:<13}{r['dg_switch']:>7.2f}{r['closed_frac']:>8.1%}"
              f"{r['kd_app_nM']:>8.0f}{r['d_kd_nM']:>+8.0f}{r['engagement']:>6.2f}"
              f"{r['selectivity']:>7.2f}{r['tether_nt']:>6}")

    if out:
        dst = HERE / "mismatches.csv"
        with dst.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0]))
            w.writeheader()
            w.writerows(out)
        print(f"\nwrote {dst} ({len(out)} rows)")
        short = min(out, key=lambda r: r["tether_nt"])
        print(f"shortest passing tether: {short['tether_nt']} nt -> {short['name']} "
              f"({short['kd_app_nM']:.0f} nM)")
    else:
        print("\nno mismatch variant passed; the filters may be too tight for this parent")


if __name__ == "__main__":
    main()
