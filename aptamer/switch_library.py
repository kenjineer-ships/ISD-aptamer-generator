"""Enumerate intramolecular-strand-displacement switches on a parent aptamer and rank them.

Construct:  5'-[aptamer]-[poly-T linker]-[displacement strand]-3'
The DS is the reverse complement of a window of the aptamer, so the grid is
(window start, DS length, linker length) -- no stochastic search needed.

Every window is scanned rather than only the binding site, because the binding site is
unknown for these aptamers and top-performing switches often hybridise at motifs outside
the known pocket (Saunders/Thompson/Soh 2025, PMC11883736).

Thermodynamics, per construct:
    G(closed) = MFE(aptamer + linker + DS)          DS paired to the aptamer
    G(open)   = MFE(aptamer) + MFE(linker + DS)     aptamer folds natively, DS does its own thing
    dG_switch = G(closed) - G(open)
    K_closed  = exp(-dG_switch / RT)
    KD_app    = KD_parent * (1 + K_closed)          switching always costs affinity
"""

import csv
import math
import pathlib
from functools import lru_cache

import primer3
import RNA
from build_parents import BLOCKS, TEMPLATE, fold, parent

# IL-6-7326.1: the truncated form the 27 nM was actually measured on -- the tightest parent
# available, and at 45 nt it keeps the finished construct in the range real E-ABs work at
# (74-mer + tether ran to 128 nt, versus 20-50 nt for working sensors). parent() guarantees
# the KD belongs to the sequence it returns. Use "IL-6-9805" for the full-length comparison.
PARENT = "IL-6-7326.1"
RT = 0.6163  # kcal/mol at 37 C
COMPLEMENT = str.maketrans("ACGT", "TGCA")

# Filter on the requirements themselves, not on a dG window: the acceptable dG depends on the
# parent's KD, so a hardcoded window silently lets constructs through when the parent changes.
KD_APP_MAX_NM = 1000.0  # affinity budget
CLOSED_MIN = 0.75       # switch must actually be closed without target
MIN_ENGAGEMENT = 0.6  # fraction of DS bases actually paired to the aptamer
SELECTIVITY_MIN = 2.0  # kcal/mol the designed duplex must beat the DS's own self-structure by
DS_LENGTHS = range(6, 15)
# Sweep out to 43 nt so the grid covers the regime the model was actually benchmarked in
# (Wilson 2019 used loops of 23-43); 3-12 alone is extrapolation. See benchmark.py.
LINKER_LENGTHS = range(3, 44)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def k_closed_from(dg_switch):
    """Closed/open equilibrium constant from dG_switch. Note the -1.

    dG_switch is G(full construct ensemble) - G(open state), and the full-construct ensemble
    contains BOTH closed and open conformations: Z_total = Z_closed + Z_open. So

        exp(-dG_switch / RT) = Z_total / Z_open = 1 + K_closed

    Taking exp(-dG/RT) as K_closed directly inflates it by exactly 1 -- a ~20% error on KD_app
    at K_closed ~ 4, and far worse for weak switches, where it also pins log10 K at 0 instead
    of letting it fall. Returns 0 when the DS does not close at all (dG_switch >= 0).
    """
    return max(math.exp(-dg_switch / RT) - 1.0, 0.0)


def randomised_positions():
    """0-indexed positions of the randomised library blocks within the 74-mer.

    Selection only varied these, so the binding site almost certainly involves them: a DS
    covering them is far likelier to gate binding than one pairing to a fixed primer region.
    """
    fixed = TEMPLATE.split("{}")
    pos, out = 0, set()
    for f, block in zip(fixed, BLOCKS):
        pos += len(f)
        out.update(range(pos, pos + block))
        pos += block
    return out


def max_homopolymer(seq):
    best = run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


@lru_cache(maxsize=None)
def open_state(tail):
    """Ensemble dG of the linker+DS tail on its own. Cached: the tail recurs across windows."""
    return fold(tail)[2]


@lru_cache(maxsize=None)
def dimer_screen(ds, window_seq):
    """(designed duplex, homodimer, hairpin) dG in kcal/mol -- primer3, 37 C, 127 mM Na+.

    The electrode carries a dense monolayer of *identical* constructs, so the effective local
    concentration of a perfectly complementary neighbour is enormous. If a DS would rather
    pair with a copy of itself than with its target window, the switch is broken before any
    IL-6 arrives. dG of 0.0 means primer3 found no such structure, which is ideal.

    Note: a BLAST screen against the genome cannot answer this -- every 6-14 nt DS occurs by
    chance tens to tens of thousands of times in 6.2 Gbp, so it flags everything and ranks
    nothing. Sequestration by a neighbour is the failure mode that actually discriminates.
    """
    kw = dict(mv_conc=127.0, dv_conc=0.0, dna_conc=50.0, temp_c=37.0)
    return (primer3.calc_heterodimer(ds, window_seq, **kw).dg / 1000,
            primer3.calc_homodimer(ds, **kw).dg / 1000,
            primer3.calc_hairpin(ds, **kw).dg / 1000)


def engagement(fc, n_aptamer, ds_start, ds_len):
    """Expected fraction of DS bases paired to the aptamer, from base-pair probabilities.

    Probabilistic rather than read off the MFE dot-bracket: a DS that is engaged in 60% of the
    ensemble and free in 40% is a real design, and binary counting scored it 1.0 or 0.0
    depending on which side of the line the single best structure fell.
    """
    bpp = fc.bpp()  # 1-indexed upper triangle; aptamer indices always precede DS indices
    expected = sum(bpp[i][j]
                   for i in range(1, n_aptamer + 1)
                   for j in range(ds_start + 1, ds_start + ds_len + 1))
    return expected / ds_len


def main():
    RNA.params_load_DNA_Mathews2004()

    aptamer, kd_parent = parent(PARENT)
    n = len(aptamer)
    g_aptamer = fold(aptamer)[2]
    randomised = randomised_positions()

    rows, scanned = [], 0
    for ds_len in DS_LENGTHS:
        for start in range(n - ds_len + 1):
            ds = revcomp(aptamer[start:start + ds_len])
            gc = (ds.count("G") + ds.count("C")) / ds_len
            if max_homopolymer(ds) > 4 or not 0.25 <= gc <= 0.75:
                continue
            for linker_len in LINKER_LENGTHS:
                tail = "T" * linker_len + ds
                construct = aptamer + tail
                scanned += 1
                ss, _, g_closed, fc = fold(construct)
                dg = g_closed - (g_aptamer + open_state(tail))
                k_closed = k_closed_from(dg)
                closed_frac = k_closed / (1 + k_closed)
                kd_app_nM = kd_parent * (1 + k_closed) * 1e9
                if closed_frac < CLOSED_MIN or kd_app_nM > KD_APP_MAX_NM:
                    continue
                eng = engagement(fc, n, n + linker_len, ds_len)
                if eng < MIN_ENGAGEMENT:
                    continue
                designed, homo, hairpin = dimer_screen(ds, aptamer[start:start + ds_len])
                selectivity = min(homo, hairpin) - designed  # >0: designed duplex wins
                rows.append(dict(
                    ds=ds, window=f"{start}-{start + ds_len - 1}", ds_len=ds_len,
                    linker_len=linker_len, gc=round(gc, 2), dg_switch=round(dg, 2),
                    closed_frac=round(closed_frac, 3), kd_app_nM=round(kd_app_nM, 1),
                    engagement=round(eng, 2),
                    rand_covered=len(randomised & set(range(start, start + ds_len))),
                    dg_designed=round(designed, 2), dg_homodimer=round(homo, 2),
                    dg_hairpin=round(hairpin, 2), selectivity=round(selectivity, 2),
                    tether_nt=linker_len + ds_len,
                    structure=ss, construct=construct,
                ))

    # Rank on things we unambiguously want, given the filters already guarantee the switch
    # closes (>=75%) and the DS is engaged (>=MIN_ENGAGEMENT):
    #   1. binding-site coverage - selection only varied the randomised positions
    #   2. sensitivity - lowest apparent KD
    #   3. speed - shortest tether
    # Engagement is deliberately NOT a ranking term. Once probabilistic it became continuous,
    # and sorting on it rewards the DS binding hardest, which is precisely what costs affinity;
    # it pushed the whole top of the table against the 1 uM ceiling.
    rows.sort(key=lambda r: (-r["rand_covered"], r["kd_app_nM"], r["tether_nt"]))
    risky = [r for r in rows if r["selectivity"] < SELECTIVITY_MIN]
    print(f"{PARENT}: {scanned} constructs scanned, {len(rows)} passed")
    print(f"dimer screen: {len(risky)}/{len(rows)} have selectivity < "
          f"{SELECTIVITY_MIN} kcal/mol (self-structure competes with the designed duplex)\n")
    print(f"{'DS':<15}{'window':>8}{'lnk':>5}{'dG':>7}{'closed':>8}"
          f"{'KDapp/nM':>10}{'eng':>6}{'rand':>6}{'selec':>7}{'tether':>8}")
    for r in rows[:15]:
        flag = " !" if r["selectivity"] < SELECTIVITY_MIN else ""
        print(f"{r['ds']:<15}{r['window']:>8}{r['linker_len']:>5}{r['dg_switch']:>7.2f}"
              f"{r['closed_frac']:>8.1%}{r['kd_app_nM']:>10.1f}"
              f"{r['engagement']:>6.2f}{r['rand_covered']:>6}{r['selectivity']:>7.2f}"
              f"{r['tether_nt']:>8}{flag}")

    fast = min(rows, key=lambda r: r["tether_nt"]) if rows else None
    if fast:
        print(f"\ntether_nt = linker + DS, an ORDINAL speed proxy (shorter = faster response).")
        print(f"Direction only, from the two published k_obs points: Wilson 2019 measured "
              f"8_23 (tether 31) at 3.58/s and 9_36 (tether 45) at 0.31/s. Two points cannot "
              f"calibrate a rate, so do not read this as seconds.")
        print(f"shortest tether among passing constructs: {fast['tether_nt']} nt "
              f"({fast['ds']}, window {fast['window']}, linker {fast['linker_len']}, "
              f"KDapp {fast['kd_app_nM']:.0f} nM)")
        span = sorted(r["tether_nt"] for r in rows)
        print(f"tether range: {span[0]}-{span[-1]} nt  "
              f"(construct length {len(aptamer) + span[0]}-{len(aptamer) + span[-1]} nt)")

    out = pathlib.Path(__file__).with_name("switches.csv")
    if rows:
        with out.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {out} ({len(rows)} rows)")
    else:
        print("\nno constructs passed; widen DG_WINDOW or lower MIN_ENGAGEMENT")


if __name__ == "__main__":
    main()
