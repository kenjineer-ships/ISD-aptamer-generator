"""Negative controls: does the scorer detect real complementarity, or just "more DNA"?

Any appended sequence lowers a construct's free energy by forming *some* structure, so a
negative dG_switch on its own proves nothing. The discriminating question is whether designed
displacement strands separate from composition-matched controls that cannot pair with the
intended window.

Three controls per designed construct, each matched on length and base composition:

  scrambled  DS bases shuffled            - same composition, complementarity destroyed
  reversed   DS reversed, not complemented - same composition, wrong pairing register
  foreign    revcomp of a window taken from a SHUFFLED aptamer - a well-formed duplex-former
             that happens not to be complementary to this aptamer

The third is the strict one: scrambled and reversed sequences may be poor duplex-formers for
generic reasons, whereas `foreign` is a genuine reverse complement of *something*. If designed
beats foreign, the score is specific to this aptamer rather than to DNA-likeness.

Reads aptamer/switches.csv, writes aptamer/negative_controls.csv
"""

import csv
import pathlib
import random
import statistics

from mismatch_tune import APTAMER, G_APTAMER, passes, score
from switch_library import revcomp

SEED = 20260815  # fixed: controls must be reproducible
N_SAMPLE = 80    # designed constructs to test, spread across the ranked file
HERE = pathlib.Path(__file__).parent
RNG = random.Random(SEED)
SHUFFLED_APTAMER = "".join(RNG.sample(APTAMER, len(APTAMER)))


def controls_for(ds, window_start):
    """One control DS of each kind, matched to `ds` on length and composition."""
    scrambled = "".join(RNG.sample(ds, len(ds)))
    reversed_ds = ds[::-1]
    start = RNG.randrange(0, len(SHUFFLED_APTAMER) - len(ds) + 1)
    foreign = revcomp(SHUFFLED_APTAMER[start:start + len(ds)])
    return {"scrambled": scrambled, "reversed": reversed_ds, "foreign": foreign}


def main():
    rows = list(csv.DictReader((HERE / "switches.csv").open()))
    step = max(1, len(rows) // N_SAMPLE)
    sample = rows[::step][:N_SAMPLE]
    print(f"{len(sample)} designed constructs sampled from {len(rows)} "
          f"(every {step}th, so the whole rank range is covered)\n")

    arms = {k: [] for k in ("designed", "scrambled", "reversed", "foreign")}
    out = []
    for base in sample:
        ds_wt = base["ds"]
        start = int(base["window"].split("-")[0])
        linker = int(base["linker_len"])

        variants = {"designed": ds_wt, **controls_for(ds_wt, start)}
        for arm, ds in variants.items():
            s = score(APTAMER, G_APTAMER, ds, start, linker)
            s.update(arm=arm, ds=ds, window=base["window"], linker_len=linker,
                     passes=bool(passes(s, ds)))
            arms[arm].append(s)
            out.append(s)

    def med(arm, key):
        return statistics.median(a[key] for a in arms[arm])

    print(f"{'arm':<11}{'n':>4}{'med dG':>9}{'med closed':>12}{'med KDapp':>11}"
          f"{'med eng':>9}{'pass rate':>11}")
    for arm in arms:
        a = arms[arm]
        rate = sum(x["passes"] for x in a) / len(a)
        print(f"{arm:<11}{len(a):>4}{med(arm, 'dg_switch'):>9.2f}"
              f"{med(arm, 'closed_frac'):>11.1%}{med(arm, 'kd_app_nM'):>11.0f}"
              f"{med(arm, 'engagement'):>9.2f}{rate:>10.0%}")

    # separation: how often does a control reach the designed median?
    d_eng = med("designed", "engagement")
    print(f"\ndesigned median engagement {d_eng:.2f}; fraction of each control reaching it:")
    for arm in ("scrambled", "reversed", "foreign"):
        frac = sum(1 for a in arms[arm] if a["engagement"] >= d_eng) / len(arms[arm])
        print(f"  {arm:<10} {frac:>5.0%}")

    verdict = all(med(arm, "engagement") < d_eng - 0.1
                  for arm in ("scrambled", "reversed", "foreign"))
    print(f"\nVERDICT: {'controls separate cleanly - score is complementarity-specific' if verdict else 'CONTROLS NOT SEPARATED - the score may be generic, investigate'}")

    dst = HERE / "negative_controls.csv"
    with dst.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"wrote {dst} ({len(out)} rows)")


if __name__ == "__main__":
    main()
