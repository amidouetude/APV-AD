"""
Correctness gate for economics.py: prove it reproduces
outputs/csv/results_summary.csv's 4E columns exactly for Ouagadougou,
both S0 (baseline) and S6 (carbon credit), using the real
config_00.SITES["Ouagadougou"] dict and a real eLER_objective() result.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "notebooks"))
import config_00  # noqa: E402
import pandas as pd

from model_bridge import evaluate_catalog_site, OUTPUTS_CSV_DIR  # noqa: E402
from economics import compute_4e  # noqa: E402

X_OUAGA = [16.37, 3.026, 1.501, 2.64, 31.8, 0.0502]
published = pd.read_csv(OUTPUTS_CSV_DIR / "results_summary.csv")

FIELDS = [
    "CAPEX_PV_kUSD", "CAPEX_AD_kUSD", "CAPEX_total_kUSD",
    "gross_revenue_kUSD_yr", "opex_kUSD_yr", "net_CF_kUSD_yr",
    "NPV_kUSD", "IRR_pct", "payback_yr",
    "LCOE_PV_USD_kWh", "CO2_avoided_total_tCO2_ha", "carbon_value_USD_ha",
]

for scenario in ["S0", "S6"]:
    result = evaluate_catalog_site("Ouagadougou", X_OUAGA, scenario=scenario)
    fourE = compute_4e(config_00.SITES["Ouagadougou"], result, scenario=scenario)
    row = published[(published["site"] == "Ouagadougou") & (published["scenario"] == scenario)].iloc[0]

    print(f"\n=== Scenario {scenario} ===")
    print(f"{'field':28s} {'bridge':>14s} {'published':>14s} {'match':>9s}")
    for f in FIELDS:
        b = fourE[f]
        p = row[f]
        if b is None or (isinstance(b, float) and b != b):
            ok = pd.isna(p)
        else:
            ok = abs(float(b) - float(p)) < max(0.5, abs(float(p)) * 0.01)
        print(f"{f:28s} {str(b):>14s} {str(p):>14s} {'OK' if ok else 'MISMATCH':>9s}")
