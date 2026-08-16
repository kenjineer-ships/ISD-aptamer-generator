"""Reconstruct the four published Neomer IL-6 aptamers and fold them (ViennaRNA, DNA params).

Sequence  = library template + per-candidate random modules (Table 2).
Kinetics  = SPR, Table 3. None bound human serum albumin.
Source: Meehan et al., PLOS ONE 2025, doi:10.1371/journal.pone.0311497
"""

import json
import pathlib

import RNA

# Library template; each {} is one randomised block. 74 nt total.
TEMPLATE = "TGTGTATAAGTC{}GAGG{}GAAT{}AACCATCGGCGCCAACA{}CATTC{}CAGA{}TCTACTAGTCAC"
BLOCKS = (2, 3, 3, 3, 3, 2)  # module A | module B
assert TEMPLATE.count("{}") == len(BLOCKS)

IL6_MW = 23718  # Da, UniProt P05231

# name -> (random modules, KD [M], koff [1/s], variant assayed).  Ordered tightest first.
APTAMERS = {
    "IL-6-7326": (("AA", "GTT", "TGC", "AAA", "AAA", "AA"), 2.70e-8, 1.50e-3, "7326.1 truncated"),
    "IL-6-6449": (("GC", "TTT", "GG", "AAA", "AAA", "AA"), 3.21e-8, 1.52e-3, "6449 full"),
    "IL-6-9805": (("AC", "TTT", "AAG", "AAA", "AAA", "AA"), 4.05e-8, 1.77e-3, "9805 full"),
    "IL-6-4202": (("CT", "TGT", "ACA", "AAA", "AAA", "AA"), 7.32e-8, 2.98e-3, "4202.1 truncated"),
}

# The paper truncated two candidates before the binding assay -- "the truncations in both
# cases effectively removed module B" -- and Table 3's KD values are for those truncated
# forms. Module B was identical (AAA/AAA/AA) across all top-15 sequences, so it carried no
# selectivity. Reconstructed as positions 0..44, i.e. through the KasI fixed block, dropping
# module B and the 3' primer site.
#
# CAVEAT: the exact truncation boundary is inferred from that sentence, not read off Fig 4
# (the FORNA panels defeated transcription). A few nt either way changes the predicted fold.
TRUNC_END = 45
TRUNCATED_OF = {"IL-6-7326.1": "IL-6-7326", "IL-6-4202.1": "IL-6-4202"}


def parent(name):
    """(sequence, KD in M) for a parent, full-length or truncated.

    The returned KD always belongs to the returned sequence. That pairing is the invariant:
    IL-6-7326's published 27 nM was measured on the truncated 7326.1, so using it with the
    full 74-mer optimises against a different molecule.
    """
    if name in TRUNCATED_OF:
        base = TRUNCATED_OF[name]
        return TEMPLATE.format(*APTAMERS[base][0])[:TRUNC_END], APTAMERS[base][1]
    mods, kd, _, variant = APTAMERS[name]
    assert "full" in variant, f"{name}'s KD was measured on '{variant}', not full-length"
    return TEMPLATE.format(*mods), kd


# Wound-fluid IL-6, pg/mL. Rembe et al. 2025, PMC11978031.
WOUND = {
    "post-surgical (median)": 986.6,
    "inflammatory (median)": 4964,
    "infected (median)": 5883,
    "max observed": 135500,
}


def _fc(seq, temp_c):
    md = RNA.md()
    md.temperature, md.noLP = temp_c, 1
    return RNA.fold_compound(seq, md)


def mfe(seq, temp_c=37.0):
    """(structure, minimum free energy). Call params_load_DNA_Mathews2004() first."""
    return _fc(seq, temp_c).mfe()


def fold(seq, temp_c=37.0):
    """(mfe_structure, mfe_dG, ensemble_dG, fold_compound).

    The ensemble free energy, not the MFE, is the correct free energy of a state: MFE reports
    only the single best structure, so a ranking built on it jumps whenever an alternative fold
    overtakes the previous winner. That discreteness was visible as non-monotonic spikes down
    the linker column. The returned fold_compound already has its partition function computed,
    so callers can take base-pair probabilities from it without refolding.
    """
    fc = _fc(seq, temp_c)
    ss, dg = fc.mfe()
    fc.exp_params_rescale(dg)  # keeps the pf numerically stable on longer constructs
    _, ens = fc.pf()
    return ss, dg, ens, fc


def occupancy(kd, pg_ml):
    """(molar concentration, fraction of aptamer bound) for an IL-6 level in pg/mL."""
    conc = pg_ml * 1e-9 / IL6_MW
    return conc, conc / (kd + conc)


def main():
    RNA.params_load_DNA_Mathews2004()

    rows = []
    for name, (mods, kd, koff, variant) in APTAMERS.items():
        seq = TEMPLATE.format(*mods)
        # Catches OCR/transcription damage in the module table: 6449's third block is short.
        bad = [f"block {i} {m!r} is {len(m)} nt, template expects {n}"
               for i, (m, n) in enumerate(zip(mods, BLOCKS)) if len(m) != n]
        ss22, dg22 = mfe(seq, 22.0)  # paper's conditions, reproduces their Fig 4
        ss37, dg37 = mfe(seq, 37.0)  # what a wound-fluid sensor actually sees
        rows.append(dict(name=name, sequence=seq, length=len(seq), KD_M=kd, koff_s=koff,
                         variant=variant, issues=bad, structure_22C=ss22, mfe_22C=dg22,
                         structure_37C=ss37, mfe_37C=dg37))

        print(f"{name}  {len(seq)} nt  KD {kd * 1e9:.1f} nM  "
              f"koff {koff:.2e}/s (1/koff {1 / koff:.0f} s)  [{variant}]")
        print(f"  {seq}")
        print(f"  22C {ss22} {dg22:+.2f}")
        print(f"  37C {ss37} {dg37:+.2f}")
        print("".join(f"  !! {b}\n" for b in bad), end="")
        print()

    best = next(iter(rows))  # APTAMERS is ordered by affinity
    print(f"Occupancy of {best['name']} (KD {best['KD_M'] * 1e9:.1f} nM) in wound fluid:")
    for label, pg in WOUND.items():
        conc, frac = occupancy(best["KD_M"], pg)
        print(f"  {label:<24} {pg:>8.0f} pg/mL = {conc * 1e12:>7.1f} pM -> {frac:6.1%} bound")

    out = pathlib.Path(__file__).with_name("parents.json")
    out.write_text(json.dumps(rows, indent=1), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
