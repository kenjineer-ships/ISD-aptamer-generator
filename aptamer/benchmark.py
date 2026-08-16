"""Benchmark the switch model against published ISD measurements.

Ground truth: Wilson, Hariri, Thompson, Eisenstein & Soh, "Independent control of the
thermodynamic and kinetic properties of aptamer switches", Nat Commun (2019), PMC6838323,
Supplementary Table 2 -- the closing equilibrium constant K_Q measured on an ATP aptamer
across displacement-strand and loop lengths, plus the fitted design law

    d(log10 K_Q) / d(L_DS) = 0.826 +/- 0.157

K_Q is the same quantity this pipeline calls K_closed = exp(-dG_switch / RT).

What this can and cannot test. Their DS *sequences* are not published in a form we could
reach, so we cannot recompute K_Q for their constructs and correlate per-construct. What we
can test is whether our model reproduces their measured *design laws* on a completely
different aptamer: the slope of log K_closed against DS length, and the sign of its dependence
on linker length. Those are the two behaviours the whole ranking rests on, so a mismatch would
invalidate the shortlist.
"""

import math
import statistics
import sys

from build_parents import fold, parent
from switch_library import PARENT, RT, k_closed_from, open_state, revcomp

# Wilson 2019 Supplementary Table 2. rows = loop length, cols = DS length -> measured K_Q.
# DS=5 is excluded: 21.4 / 0.0379 / 9.2 is non-monotonic in both directions, i.e. K_Q too
# small to fit rather than a real measurement.
MEASURED_KQ = {
    23: {9: 39.0, 8: 3.58, 7: 0.382, 6: 0.0314},
    25: {9: 12.8, 8: 2.29, 7: 0.272, 6: 0.0265},
    33: {9: 4.99, 8: 0.735, 7: 0.0985, 6: 0.0182},
}
PUBLISHED_SLOPE = (0.826, 0.157)  # d(log10 K_Q)/d(L_DS), value and stated error

DS_LENGTHS = range(6, 13)
LINKER_LENGTHS = (23, 25, 33)  # match theirs so the comparison is like-for-like


def fit_slope(xs, ys):
    """Least-squares slope of ys on xs."""
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    return num / den


def spearman(a, b):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos + 1.0
        return r

    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = statistics.fmean(ra), statistics.fmean(rb)
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = math.sqrt(sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb))
    return num / den if den else float("nan")


ENGAGED_MIN_LOGK = 0.1  # log10 K_closed below this = the DS never really closes


def predicted_log_kq(aptamer, g_aptamer, ds_len, linker_len, engaged_only=False):
    """median log10 K_closed over DS windows of this length.

    Median, not mean: window choice is arbitrary here (their binding site is known, ours is
    not) and a few windows collide with the aptamer's own hairpin, producing outliers.

    engaged_only drops windows where the DS does not measurably close. Wilson's DS targets the
    binding site by construction and always engages, so non-engaging windows are not a
    like-for-like comparison -- on a short parent they are numerous and pin the median at 0,
    flattening the slope. Reported alongside the unfiltered number, never instead of it.
    """
    vals = []
    for start in range(len(aptamer) - ds_len + 1):
        tail = "T" * linker_len + revcomp(aptamer[start:start + ds_len])
        # ensemble dG throughout - mixing MFE for the closed state with an ensemble open state
        # would compare different quantities
        dg = fold(aptamer + tail)[2] - (g_aptamer + open_state(tail))
        # log10 of the true K_closed. -dg/RT/ln10 would be log10(1 + K_closed), which saturates
        # at 0 for weak switches and artificially flattens the slope this test measures.
        k = k_closed_from(dg)
        vals.append(math.log10(k) if k > 0 else float("-inf"))
    if engaged_only:
        kept = [v for v in vals if v >= ENGAGED_MIN_LOGK]
        if kept:
            vals = kept
    finite = [v for v in vals if math.isfinite(v)]
    return statistics.median(finite) if finite else float("nan")


def main():
    import RNA
    RNA.params_load_DNA_Mathews2004()

    # --parent NAME compares a different parent than the pipeline's current one; the whole
    # open question here is how much of the agreement is parent-specific.
    name = PARENT
    if "--parent" in sys.argv:
        name = sys.argv[sys.argv.index("--parent") + 1]
    aptamer, _ = parent(name)
    g_aptamer = fold(aptamer)[2]
    globals()["PARENT_USED"] = name

    # Default is ALL windows. --engaged-only was tried as a fix for the truncated parent's
    # poor agreement and made it worse (slope -0.166, Spearman -0.448): conditioning on
    # "windows where the DS wins" is selection bias, since more windows qualify as DS grows,
    # so the median's composition shifts with the very variable being regressed on.
    engaged = "--engaged-only" in sys.argv
    grid = {ln: {d: predicted_log_kq(aptamer, g_aptamer, d, ln, engaged_only=engaged)
                 for d in DS_LENGTHS} for ln in LINKER_LENGTHS}
    print(f"aggregation: {'engaging windows only' if engaged else 'ALL windows'}"
          f"  (toggle with --all-windows)")

    print(f"parent {name} ({len(aptamer)} nt); predicted log10 K_closed, "
          f"median over all DS windows\n")
    hdr = "".join(f"{d:>8}" for d in DS_LENGTHS)
    print(f"{'loop':>5}{hdr}")
    for ln in LINKER_LENGTHS:
        print(f"{ln:>5}" + "".join(f"{grid[ln][d]:>8.2f}" for d in DS_LENGTHS))

    # --- test 1: slope of log K vs DS length, against their fitted design law ---
    print("\nTEST 1  d(log10 K)/d(L_DS)   published 0.826 +/- 0.157")
    lo, hi = PUBLISHED_SLOPE[0] - PUBLISHED_SLOPE[1], PUBLISHED_SLOPE[0] + PUBLISHED_SLOPE[1]
    ours = []
    for ln in LINKER_LENGTHS:
        ds = [d for d in DS_LENGTHS if d in MEASURED_KQ[ln]] or list(DS_LENGTHS)
        s_pred = fit_slope(ds, [grid[ln][d] for d in ds])
        s_meas = fit_slope(ds, [math.log10(MEASURED_KQ[ln][d]) for d in ds])
        ours.append(s_pred)
        mark = "in range" if lo <= s_pred <= hi else "OUT OF RANGE"
        print(f"  loop {ln:>2}: predicted {s_pred:+.3f}   measured {s_meas:+.3f}   [{mark}]")
    print(f"  mean predicted slope {statistics.fmean(ours):+.3f}")

    # --- test 2: sign of the linker-length dependence ---
    print("\nTEST 2  longer loop must lower K (weaker intramolecular competition)")
    for d in DS_LENGTHS:
        s_pred = fit_slope(list(LINKER_LENGTHS), [grid[ln][d] for ln in LINKER_LENGTHS])
        meas = [MEASURED_KQ[ln][d] for ln in LINKER_LENGTHS if d in MEASURED_KQ[ln]]
        s_meas = (fit_slope(list(LINKER_LENGTHS), [math.log10(v) for v in meas])
                  if len(meas) == len(LINKER_LENGTHS) else float("nan"))
        ok = "ok" if s_pred < 0 else "WRONG SIGN"
        print(f"  DS {d:>2}: predicted {s_pred:+.4f}/nt   measured {s_meas:+.4f}/nt   [{ok}]")

    # --- test 3: rank agreement across the shared grid ---
    pairs = [(ln, d) for ln in LINKER_LENGTHS for d in DS_LENGTHS if d in MEASURED_KQ[ln]]
    pred = [grid[ln][d] for ln, d in pairs]
    meas = [math.log10(MEASURED_KQ[ln][d]) for ln, d in pairs]
    print(f"\nTEST 3  Spearman(predicted, measured) over {len(pairs)} shared grid points: "
          f"{spearman(pred, meas):+.3f}")
    print("        (different aptamer and unknown DS windows, so rank agreement is the "
          "meaningful statistic, not absolute values)")


if __name__ == "__main__":
    main()
