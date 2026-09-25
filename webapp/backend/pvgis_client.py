"""
Live PVGIS-SARAH3 TMY data ingestion for arbitrary points, producing the
same DataFrame shape the real (unmodified) simulation_functions.py
pipeline expects -- mirroring whatever notebooks/01_load_data.ipynb does
for the four catalog sites (see model_bridge.load_hourly_csv's docstring
for the DatetimeIndex requirement this must also satisfy).
"""
import urllib.request
import json
import pandas as pd

PVGIS_TMY_URL = "https://re.jrc.ec.europa.eu/api/v5_3/tmy"

# PVGIS's own field names -> the names simulation_functions.py's pipeline
# expects (matches outputs/csv/hourly_{site}.csv exactly).
FIELD_MAP = {
    "G(h)": "GHI",
    "Gb(n)": "DNI",
    "Gd(h)": "DHI",
    "IR(h)": "IR",
    "WS10m": "WS",
    "WD10m": "WD",
}


def fetch_tmy(lat: float, lon: float, timeout: int = 30) -> pd.DataFrame:
    """
    Fetch a PVGIS-SARAH3 TMY for (lat, lon). Returns an hourly DataFrame
    with a real DatetimeIndex plus T2m, RH, GHI, DNI, DHI, IR, WS, WD, SP,
    month, hour_of_day, day_of_year. Does NOT set in_season -- see
    add_growing_season() below, a growing-season assumption is a separate,
    explicit step, never silently defaulted here.

    PVGIS splices each calendar month from a different "typical" year
    (e.g. Jan from 2005, Jun from 2022) -- each date's own day-of-year is
    still correct for solar-position purposes to within ~1 day around a
    leap-year boundary, which is not corrected for here.
    """
    url = f"{PVGIS_TMY_URL}?lat={lat}&lon={lon}&outputformat=json"
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        payload = json.loads(resp.read())

    meta = payload["inputs"]["meteo_data"]
    rows = payload["outputs"]["tmy_hourly"]
    df = pd.DataFrame(rows)
    df = df.rename(columns=FIELD_MAP)
    df["time_UTC"] = pd.to_datetime(df["time(UTC)"], format="%Y%m%d:%H%M")
    df = df.drop(columns=["time(UTC)"]).set_index("time_UTC").sort_index()

    df["month"] = df.index.month
    df["hour_of_day"] = df.index.hour
    df["day_of_year"] = df.index.dayofyear

    df.attrs["radiation_db"] = meta["radiation_db"]
    df.attrs["elevation_m"] = payload["inputs"]["location"]["elevation"]
    return df


def add_growing_season(df: pd.DataFrame, planting_mmdd: str, harvest_mmdd: str,
                        source: str = "assumed default -- not a sourced crop calendar") -> pd.DataFrame:
    """
    Flag in_season = True for rows within [planting_mmdd, harvest_mmdd]
    ("MM-DD", inclusive; handles a season that wraps across Dec 31/Jan 1).

    This is the crude fallback from docs/webapp_africa_platform_design.md
    section 2 -- NOT the FAO/GIEWS crop-calendar lookup described there.
    `source` is always stamped into df.attrs so a caller can never lose
    track of which method produced a given in_season flag.
    """
    df = df.copy()
    ref_year = 2001  # non-leap reference year, only used for day-of-year math
    p_doy = pd.Timestamp(f"{ref_year}-{planting_mmdd}").dayofyear
    h_doy = pd.Timestamp(f"{ref_year}-{harvest_mmdd}").dayofyear
    if p_doy <= h_doy:
        in_season = (df["day_of_year"] >= p_doy) & (df["day_of_year"] <= h_doy)
    else:
        in_season = (df["day_of_year"] >= p_doy) | (df["day_of_year"] <= h_doy)
    df["in_season"] = in_season.values
    df.attrs["season_source"] = source
    df.attrs["planting_mmdd"] = planting_mmdd
    df.attrs["harvest_mmdd"] = harvest_mmdd
    return df
