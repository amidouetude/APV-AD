"""
Correctness gate for pareto.py: run the real NSGA-II search (via the
runtime-injection path, throwaway key -- NOT the real "Konya" catalog
entry) against Konya's own real hourly data, and compare the resulting
front's extremes to the manuscript's published values (Section 5.3,
"Beyond the scalar optimum"): water-maximising end (eLER=1.20,
LER_water=36.1%) to eLER-maximising end (eLER=1.64, LER_water=18.9%).
"""
import copy
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "notebooks"))
import config_00  # noqa: E402

from model_bridge import load_hourly_csv  # noqa: E402
from pareto import run_pareto  # noqa: E402

site_dict = copy.deepcopy(config_00.SITES["Konya"])
df = load_hourly_csv("Konya")

t0 = time.time()


def progress(gen, n_gen):
    if gen % 10 == 0:
        print(f"  gen {gen:3d}/{n_gen}  elapsed = {time.time()-t0:.0f}s")


eLER, LER_water, X = run_pareto("_validate_konya", site_dict, df,
                                 pop_size=50, n_gen=80, progress_cb=progress)
elapsed = time.time() - t0

order = eLER.argsort()
eLER_s, water_s = eLER[order], LER_water[order]

print()
print(f"Runtime: {elapsed:.0f}s ({len(eLER)} front solutions)")
print()
print(f"eLER range on front:      {eLER.min():.3f} -- {eLER.max():.3f}   (published: 1.20 -- 1.64)")
print(f"LER_water range on front: {water_s.min():.1f}% -- {water_s.max():.1f}%   (published: 18.9 -- 36.1)")
print()
print("Water-maximising end (lowest eLER):")
print(f"  eLER = {eLER_s[0]:.3f}  LER_water = {water_s[0]:.1f}%   (published: eLER=1.20, LER_water=36.1%)")
print("eLER-maximising end (highest eLER):")
print(f"  eLER = {eLER_s[-1]:.3f}  LER_water = {water_s[-1]:.1f}%   (published: eLER=1.64, LER_water=18.9%)")
