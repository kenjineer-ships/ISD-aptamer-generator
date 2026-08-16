import sys, time, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_cofold as rc

t = rc.fetch_targets()
d = rc.load_dna_chains()
t0 = time.time()
rec = rc.run_one("boltz2", d[0], t["IL6"])
print("WALL", round(time.time() - t0, 1))
print("metrics", json.dumps(rec["metrics"], indent=2)[:2000])
print("epitope n", rec.get("epitope", {}).get("n_contact_residues"))
print("chains", rec.get("chain_ids_in_output"))
(rc.RESULTS / (rec["key"] + ".json")).write_text(json.dumps(rec, indent=2))
print("saved")
