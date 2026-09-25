"""
Non-invasive bridge to the validated APV+AD simulation engine.

This module does NOT modify notebooks/simulation_functions.py or
notebooks/config_00.py in any way -- it only imports them and makes
arbitrary (non-catalog) locations usable by the same unmodified
eLER_objective()/compute_AD(), which look up site parameters via
SITES[site].

CONCURRENCY (corrected 2026-09)
-------------------------------
The previous mechanism wrote the runtime site straight into the shared
module-level dict config_00.SITES and popped it in a finally block. The
phase-0 audit (audit/phase0/) measured what that does under load, using
the exact key api.py builds -- f"_runtime_{lat}_{lon}", which depends on
the coordinates ALONE while the site dict also depends on the crop and on
borrow_from:

    same coordinates, different crops   : 12 threads -> 11 KeyError, 1 silently
                                          contaminated result, 0 correct
    same coordinates, same parameters   : 10 threads -> 9 KeyError, 1 correct
    different coordinates (control)     :  8 threads -> 8 correct

Two distinct failures: one thread's finally-pop removes the key while
another is still mid-pipeline (KeyError -> HTTP 500), and two requests on
the same point with different parameters overwrite each other's entry
(silently wrong result). Note the control: a concurrency test that varies
the coordinates passes and detects neither.

The fix keeps the shared dict read-only and gives every request its own
overlay, scoped by contextvars. Threads do not inherit a context, and each
asyncio task carries its own, so FastAPI isolates requests either way --
whether the endpoint is sync (thread-pool) or async (event loop).
"""
import contextvars
import sys
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]          # .../APV_AD_Project
NOTEBOOKS_DIR = ROOT / "notebooks"
OUTPUTS_CSV_DIR = ROOT / "outputs" / "csv"

sys.path.insert(0, str(NOTEBOOKS_DIR))

import config_00                       # noqa: E402  (unmodified, real config)
import simulation_functions as sim     # noqa: E402  (unmodified, real pipeline)

import pandas as pd                    # noqa: E402


# ----------------------------------------------------------------------
# Per-request site overlay
# ----------------------------------------------------------------------
_RUNTIME_SITES: contextvars.ContextVar = contextvars.ContextVar(
    "apvad_runtime_sites", default=None
)


class UnknownSiteError(KeyError):
    """
    Raised when a site key resolves neither in the catalog nor in the
    current request's overlay.

    Subclasses KeyError so any existing `except KeyError` still catches it,
    but carries a message an API layer can surface as a 4xx instead of
    letting a bare KeyError escape as an undocumented HTTP 500.
    """


class _SiteRegistry(dict):
    """
    The four catalog sites, plus a per-context overlay.

    Reads consult the overlay first, then the shared catalog. Nothing ever
    writes to the shared catalog at request time, so two concurrent requests
    cannot see or clobber each other's entry.

    keys()/items()/iteration deliberately expose the CATALOG ONLY: callers
    such as api.CATALOG_SITES enumerate this to list the published reference
    sites, and a throwaway runtime entry has no business appearing there.
    """

    def _overlay(self):
        return _RUNTIME_SITES.get()

    def __getitem__(self, key):
        overlay = self._overlay()
        if overlay is not None and key in overlay:
            return overlay[key]
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            raise UnknownSiteError(
                f"Unknown site '{key}'. Known catalog sites: "
                f"{sorted(dict.keys(self))}. A runtime site is only visible "
                f"inside the runtime_site() context that registered it."
            ) from None

    def __contains__(self, key):
        overlay = self._overlay()
        if overlay is not None and key in overlay:
            return True
        return dict.__contains__(self, key)

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default


def install_registry() -> _SiteRegistry:
    """
    Swap config_00.SITES for the context-aware registry, once.

    simulation_functions.py does `from config_00 import SITES, ...` at import
    time, so it holds its OWN reference to the original dict. Rebinding
    config_00.SITES alone would leave the pipeline reading the old object,
    and the overlay would be silently ignored -- so sim.SITES is rebound to
    the same registry object here. Both names must point at it.
    """
    if isinstance(config_00.SITES, _SiteRegistry):
        return config_00.SITES
    registry = _SiteRegistry(config_00.SITES)
    config_00.SITES = registry
    sim.SITES = registry
    return registry


SITES = install_registry()


@contextmanager
def runtime_site(key: str, site_dict: dict):
    """
    Make `key` resolve to `site_dict` for the duration of this block, in
    THIS context only.

    Replaces the old register_runtime_site()/pop() pair. Nested use is
    supported: the overlay is copied, not mutated in place, so an inner
    block cannot corrupt an outer one.
    """
    current = _RUNTIME_SITES.get()
    overlay = dict(current) if current else {}
    overlay[key] = site_dict
    token = _RUNTIME_SITES.set(overlay)
    try:
        yield key
    finally:
        _RUNTIME_SITES.reset(token)


def register_runtime_site(key: str, site_dict: dict) -> None:
    """
    DEPRECATED -- kept so existing callers keep working, but it now writes
    to the per-context overlay rather than the shared catalog, and the entry
    lives until the surrounding context ends rather than until someone pops
    it. Prefer the runtime_site() context manager, which has an explicit
    scope.
    """
    current = _RUNTIME_SITES.get()
    overlay = dict(current) if current else {}
    overlay[key] = site_dict
    _RUNTIME_SITES.set(overlay)


# ----------------------------------------------------------------------
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


def evaluate_point(key: str, site_dict: dict, hourly_df: pd.DataFrame, x,
                    scenario: str = "S0", **kwargs) -> dict:
    """
    Evaluate the real pipeline at an arbitrary point, with the site
    definition visible only to this request.

    `site_dict` must carry every key SITES[...] normally provides (lat, lon,
    min_beta_deg, crop/substrate/economic parameters -- see config_00.SITES
    for the full shape); `hourly_df` must have the same columns as the
    catalog hourly_{site}.csv files (GHI, DNI, DHI, T2m, WS, in_season, ...).

    `lat` may be negative: since the 2026-09 hemisphere correction the
    pipeline orients the collector towards the equator on its own.
    """
    with runtime_site(key, site_dict):
        return sim.eLER_objective(x, hourly_df, site=key, scenario=scenario,
                                   return_full=True, **kwargs)
