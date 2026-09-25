"""
Correctness gate (see docs/webapp_africa_platform_design.md, section 4):
prove the bridge reproduces the manuscript's published Ouagadougou / S0
row (outputs/csv/results_summary.csv) exactly, using the real, unmodified
simulation_functions.py + config_00.py -- no shortcuts, no re-derivation.
"""
import pandas as pd
from model_bridge import evaluate_catalog_site, OUTPUTS_CSV_DIR

# Optimized design vector for Ouagadougou, S0, taken verbatim from
# outputs/csv/results_summary.csv (beta, d_row, H_m, V_dig, HRT, f_PV_heat).
X_OUAGA = [16.37, 3.026, 1.501, 2.64, 31.8, 0.0502]

result = evaluate_catalog_site("Ouagadougou", X_OUAGA, scenario="S0")

published = pd.read_csv(OUTPUTS_CSV_DIR / "results_summary.csv")
row = published[(published["site"] == "Ouagadougou") & (published["scenario"] == "S0")].iloc[0]

CHECKS = [
    "LER_crop", "LER_PV_defA", "LER_biogas", "eLER_defA", "LER_2C",
    "PV_sold_MWh_ha", "biogas_total_MWh_ha", "ESR",
    "GCR_pct", "OLR_kgVS_m3d",
    "W_saved_mm_season", "LER_water_season_pct",
]

print(f"{'field':22s} {'bridge':>12s} {'published':>12s} {'match':>9s}")
all_ok = True
for key in CHECKS:
    b = float(result[key])
    p = float(row[key])
    ok = abs(b - p) < 0.02
    all_ok &= ok
    print(f"{key:22s} {b:12.4f} {p:12.4f} {'OK' if ok else 'MISMATCH':>9s}")

print()
print("RESULT: ALL MATCH" if all_ok else "RESULT: MISMATCH FOUND")
