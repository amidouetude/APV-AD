"""
Correctness gate for optimizer.py: run the real DE + L-BFGS-B optimizer
(via the runtime-injection path, using a throwaway key -- NOT the real
"Ouagadougou" catalog entry) against Ouagadougou's own real hourly data,
and compare the result to the manuscript's published optimum
(outputs/csv/results_summary.csv, Ouagadougou/S0) as a check that the
optimizer wiring itself -- not just single-point evaluation -- is correct.
"""
import copy
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "notebooks"))
import config_00  # noqa: E402
import pandas as pd

from model_bridge import load_hourly_csv, OUTPUTS_CSV_DIR  # noqa: E402
from optimizer import optimize_site  # noqa: E402

site_dict = copy.deepcopy(config_00.SITES["Ouagadougou"])
df = load_hourly_csv("Ouagadougou")

t0 = time.time()


def progress(gen, best):
    if gen % 10 == 0:
        print(f"  gen {gen:3d}/50  best eLER so far = {best:.4f}  "
              f"elapsed = {time.time()-t0:.0f}s")


x_opt, full, history = optimize_site("_validate_ouaga", site_dict, df, scenario="S0",
                                      progress_cb=progress)
elapsed = time.time() - t0

published = pd.read_csv(OUTPUTS_CSV_DIR / "results_summary.csv")
row = published[(published["site"] == "Ouagadougou") & (published["scenario"] == "S0")].iloc[0]

print()
print(f"Runtime: {elapsed:.0f}s  (paper reports 102-172s per site)")
print()
print(f"{'field':16s} {'optimizer found':>16s} {'published':>12s}")
pairs = [
    ("beta_deg", full["beta_deg"], row["beta_deg"]),
    ("d_row_m", full["d_row_m"], row["d_row_m"]),
    ("H_m_m", full["H_m_m"], row["H_m_m"]),
    ("V_dig_m3", full["V_dig_m3"], row["V_dig_m3"]),
    ("HRT_days", full["HRT_days"], row["HRT_days"]),
    ("f_PV_heat", full["f_PV_heat"], row["f_PV_heat"]),
    ("eLER_defA", full["eLER_defA"], row["eLER_defA"]),
]
for name, mine, pub in pairs:
    print(f"{name:16s} {mine:16.4f} {pub:12.4f}")

print()
print(f"Convergence: gen 1 eLER = {history[0]:.4f}  ->  gen {len(history)} eLER = {history[-1]:.4f}")
print()
print("Boundary-convergence check (paper Finding 1: d_row -> ~3.03m, H_m -> ~1.50m):")
print(f"  d_row* = {full['d_row_m']:.3f} m  (published minimum: 3.026 m)")
print(f"  H_m*   = {full['H_m_m']:.3f} m  (published minimum: 1.501 m)")
