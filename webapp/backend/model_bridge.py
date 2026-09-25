"""
Non-invasive bridge to the validated APV+AD simulation engine.

This module does NOT modify notebooks/simulation_functions.py or
notebooks/config_00.py in any way -- it only imports them and, for
arbitrary (non-catalog) locations, registers a new entry into the
in-memory config_00.SITES dict at call time, never touching either
file on disk. This is the mechanism for "any African point" described
in docs/webapp_africa_platform_design.md (section 2), chosen instead
of the site: str | dict signature change proposed there, specifically
to honor "don't touch the existing files."
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]          # .../APV_AD_Project
NOTEBOOKS_DIR = ROOT / "notebooks"
OUTPUTS_CSV_DIR = ROOT / "outputs" / "csv"

sys.path.insert(0, str(NOTEBOOKS_DIR))

import config_00                       # noqa: E402  (unmodified, real config)
import simulation_functions as sim     # noqa: E402  (unmodified, real pipeline)

import pandas as pd


def load_hourly_csv(site_name: str) -> pd.DataFrame:
    """
    Load the pre-processed hourly climate CSV for one of the 4 catalog sites.

    compute_solar_angles() requires a real pd.DatetimeIndex (see its
    docstring: "df : hourly DataFrame with DatetimeIndex (UTC)") -- it reads
    idx.day_of_year / idx.hour off df.index, NOT off the hour_of_day/
    day_of_year columns already present in the CSV. A plain pd.read_csv()
    leaves the default RangeIndex, which pd.DatetimeIndex() silently
    reinterprets as nanoseconds-since-epoch -- every row collapses to the
    same instant (1970-01-01), so solar elevation comes out constant and
    F_shad is 0 for every hour. This mirrors whatever 01_load_data.ipynb
    does before handing off to the rest of the pipeline; it is a loader
    responsibility, not something simulation_functions.py enforces itself.
    """
    path = OUTPUTS_CSV_DIR / f"hourly_{site_name}.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"No processed hourly CSV for '{site_name}' at {path}. "
            "Run 01_load_data.ipynb for this site first."
        )
    df = pd.read_csv(path)
    df["time_UTC"] = pd.to_datetime(df["time_UTC"])
    df = df.set_index("time_UTC")
    return df


def evaluate_catalog_site(site_name: str, x, scenario: str = "S0", **kwargs) -> dict:
    """
    Evaluate the real, unmodified pipeline for one of the 4 sites already
    registered in config_00.SITES (Konya, Almeria, Ouagadougou, Freiburg).
    This is the correctness gate: it must reproduce
    outputs/csv/results_summary.csv exactly for a given (site, scenario).
    """
    df = load_hourly_csv(site_name)
    return sim.eLER_objective(x, df, site=site_name, scenario=scenario,
                               return_full=True, **kwargs)


def register_runtime_site(key: str, site_dict: dict) -> None:
    """
    Inject a new site definition into config_00.SITES for the lifetime of
    this process, WITHOUT writing to config_00.py on disk. This is how an
    arbitrary African point (not one of the 4 published sites) becomes
    usable by the exact same, unmodified eLER_objective()/compute_AD()
    functions, which look up site parameters via SITES[site].

    NOTE: this mutates a module-level dict shared across the whole process.
    Safe for a single-request-at-a-time script or dev server; a real
    multi-worker API must give each request its own SITES-equivalent
    (e.g. a per-request dict passed explicitly) rather than sharing this
    global -- flagged here rather than solved, since it's a concurrency
    concern for the eventual FastAPI service, not for this validation step.
    """
    config_00.SITES[key] = site_dict


def evaluate_point(key: str, site_dict: dict, hourly_df: pd.DataFrame, x,
                    scenario: str = "S0", **kwargs) -> dict:
    """
    Evaluate the real pipeline at an arbitrary point by registering a
    throwaway site entry, then calling the same unmodified eLER_objective().
    `site_dict` must carry every key SITES[...] normally provides (lat, lon,
    min_beta_deg, crop/substrate/economic parameters -- see config_00.SITES
    for the full shape); `hourly_df` must have the same columns as the
    catalog hourly_{site}.csv files (GHI, DNI, DHI, T2m, WS, in_season, ...).
    """
    register_runtime_site(key, site_dict)
    try:
        return sim.eLER_objective(x, hourly_df, site=key, scenario=scenario,
                                   return_full=True, **kwargs)
    finally:
        config_00.SITES.pop(key, None)
