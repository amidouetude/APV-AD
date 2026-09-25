"""
Proof-of-concept: run the real, unmodified APV+AD pipeline on a genuinely
new African location -- Bamako, Mali -- using LIVE PVGIS-SARAH3 data
fetched over the network, not one of the four sites baked into
config_00.SITES.

Honesty rule (matches the confidence-tier system in
webapp/design/Analyze.dc.html and docs/webapp_africa_platform_design.md
section 3): only the climate/geometry side is genuinely site-specific
here (live TMY, real lat/lon/elevation, the real optimized-design-vector
physics). The crop/substrate/economic assumptions are BORROWED from
Ouagadougou (nearest published site, similar Sahelian climate) because
Bamako has no sourced crop calendar, BMP data, or country economics yet
-- this script prints that distinction explicitly rather than blurring it.
"""
import copy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "notebooks"))
import config_00  # noqa: E402

from model_bridge import evaluate_point  # noqa: E402
from pvgis_client import fetch_tmy, add_growing_season  # noqa: E402

BAMAKO_LAT, BAMAKO_LON = 12.65, -8.00

print(f"Fetching live PVGIS-SARAH3 TMY for Bamako, Mali ({BAMAKO_LAT}N, {BAMAKO_LON}E)...")
df = fetch_tmy(BAMAKO_LAT, BAMAKO_LON)
print(f"  -> {len(df)} hourly rows, radiation DB = {df.attrs['radiation_db']}, "
      f"elevation = {df.attrs['elevation_m']} m")

df = add_growing_season(df, planting_mmdd="07-01", harvest_mmdd="10-30",
                         source="BORROWED from Ouagadougou's calendar -- not a Mali-specific source")

site_dict = copy.deepcopy(config_00.SITES["Ouagadougou"])
site_dict["lat"] = BAMAKO_LAT
site_dict["lon"] = BAMAKO_LON
site_dict["elevation_m"] = df.attrs["elevation_m"]
# Everything else in site_dict (crop, substrate, economics, min_beta_deg,
# koppen, timezone) is still Ouagadougou's -- borrowed, not sourced for
# Bamako. Flagged here, not silently reused.
BORROWED_FIELDS = [k for k in site_dict if k not in ("lat", "lon", "elevation_m")]

X_STARTING_POINT = [16.37, 3.026, 1.501, 2.64, 31.8, 0.0502]  # Ouagadougou's optimized design, reused as-is

result = evaluate_point("Bamako_demo", site_dict, df, X_STARTING_POINT, scenario="S0")

print()
print("=" * 70)
print("Bamako, Mali -- real live-TMY run through the unmodified pipeline")
print("=" * 70)
print(f"{'LER_crop':22s} {result['LER_crop']:.3f}   [genuinely site-specific: live TMY + geometry]")
print(f"{'LER_PV (Def A)':22s} {result['LER_PV_defA']:.3f}   [genuinely site-specific]")
print(f"{'LER_biogas':22s} {result['LER_biogas']:.3f}   [borrowed crop/substrate assumptions]")
print(f"{'eLER (Def A)':22s} {result['eLER_defA']:.3f}")
print(f"{'PV sold (MWh/ha/yr)':22s} {result['PV_sold_MWh_ha']:.1f}  [genuinely site-specific]")
print(f"{'Biogas (MWh/ha/yr)':22s} {result['biogas_total_MWh_ha']:.2f}  [borrowed BMP/manure assumptions]")
print(f"{'GCR':22s} {result['GCR_pct']:.1f}%")
print(f"{'W_saved (mm/season)':22s} {result['W_saved_mm_season']:.1f}  [borrowed growing-season window]")
print(f"{'LER_water (season %)':22s} {result['LER_water_season_pct']:.1f}%")
print()
print("Design vector used: Ouagadougou's *optimized* geometry, reused as a starting")
print("point for Bamako -- NOT re-optimized for this site (no DE run yet).")
print()
print("Borrowed (not Bamako-specific) input fields:")
print(" ", ", ".join(BORROWED_FIELDS))
