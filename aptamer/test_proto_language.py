"""proto-language MCMC + a ViennaRNA-backed custom constraint, on Windows.

proto-tools can't dispatch tools on Windows (micromamba), but the language layer is pure
Python, so a constraint calling ViennaRNA directly never touches it. This proves the whole
switch-design loop runs natively. Exits non-zero on failure.
"""

import RNA
from build_parents import APTAMERS, TEMPLATE, mfe
from proto_language.core import Constraint, ConstraintOutput, Construct, Program, Segment
from proto_language.generator import RandomNucleotideGenerator, RandomNucleotideGeneratorConfig
from proto_language.optimizer import MCMCOptimizer, MCMCOptimizerConfig
from proto_tools.transforms.masking import MaskingStrategy

RNA.params_load_DNA_Mathews2004()

APTAMER = TEMPLATE.format(*APTAMERS["IL-6-7326"][0])
APO = mfe(APTAMER)[1]
TARGET_DG = -7.5  # kcal/mol of displacement-strand stabilisation


def switch_energy(input_sequences, config):
    """Score |dG_switch - target|, where dG_switch = MFE(aptamer+linker+DS) - MFE(aptamer)."""
    out = []
    for seqs in input_sequences:
        dg = mfe("".join(s.sequence for s in seqs))[1] - APO
        out.append(ConstraintOutput(score=min(abs(dg - TARGET_DG) / 5, 1.0),
                                    metadata={"dG_switch": round(dg, 2)}))
    return out


aptamer = Segment(sequence=APTAMER, sequence_type="dna", label="aptamer")
linker = Segment(length=6, sequence_type="dna", label="linker")
ds = Segment(length=10, sequence_type="dna", label="ds")

# One generator per independently-varying segment: assign() REPLACES its segment list, and
# passing several at once *ties* them (shared values, equal length required).
cfg = RandomNucleotideGeneratorConfig(masking_strategy=MaskingStrategy(num_mutations=2))
gens = []
for seg in (linker, ds):
    g = RandomNucleotideGenerator(cfg)
    g.assign(seg)
    gens.append(g)

construct = Construct(segments=[aptamer, linker, ds])
opt = MCMCOptimizer(
    constructs=[construct],
    generators=gens,
    constraints=[Constraint(inputs=[aptamer, linker, ds], function=switch_energy,
                            function_config={}, weight=1.0, label="switch_dG")],
    config=MCMCOptimizerConfig(num_steps=40, num_results=3, proposals_per_result=4),
)
Program(optimizers=[opt], num_results=3).run()

print(f"\napo MFE {APO:+.2f}, target dG_switch {TARGET_DG:+.1f}\n")
best = []
for i, seq in enumerate(construct.joined_sequences, 1):
    dg = mfe(seq.sequence)[1] - APO
    best.append(dg)
    print(f"{i}. linker+DS = {seq.sequence[len(APTAMER):]}   dG_switch {dg:+.2f}")

assert best, "optimizer returned no results"
assert min(abs(dg - TARGET_DG) for dg in best) < 2.0, f"MCMC did not approach target: {best}"
print("\nPROTO-LANGUAGE ON WINDOWS: OK")
