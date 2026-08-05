import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Ensure the notebooks/ folder (where simulation_functions.py is exported) is importable
ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = ROOT / "notebooks"
sys.path.insert(0, str(NOTEBOOKS))

from simulation_functions import (
    compute_solar_angles,
    compute_shading,
    compute_POA,
    modules_per_ha,
)


def make_hourly_df():
    idx = pd.date_range("2020-06-01", periods=24, freq="H")
    df = pd.DataFrame(index=idx)
    # Simple synthetic irradiance pattern (0 at night, peak midday)
    hours = np.arange(24)
    GHI = np.maximum(0, 1000 * np.sin(np.pi * (hours / 24)))
    DNI = GHI * 0.7
    DHI = GHI - DNI
    df["GHI"] = GHI
    df["DNI"] = DNI
    df["DHI"] = DHI
    df["T2m"] = 20.0
    df["WS"] = 2.0
    return df


def test_compute_solar_angles_basic():
    df = make_hourly_df()
    out = compute_solar_angles(df.copy(), lat_deg=45.0)
    assert "solar_elevation_deg" in out.columns
    vals = out["solar_elevation_deg"].values
    assert np.all(np.isfinite(vals))
    assert vals.min() >= -90.0 and vals.max() <= 90.0


def test_compute_shading_and_gcr():
    df = make_hourly_df()
    df = compute_solar_angles(df, lat_deg=0.0)
    out = compute_shading(df, beta_deg=20.0, d_row_m=4.0, H_m_m=2.0)
    assert "F_shad" in out.columns and "GCR" in out.columns
    # shading fraction should be within [0, 1]
    f = out["F_shad"].values
    assert np.all(f >= 0.0) and np.all(f <= 1.0)
    # GCR should be a finite scalar repeated per-row
    gcr_vals = out["GCR"].unique()
    assert len(gcr_vals) == 1
    assert np.isfinite(gcr_vals[0]) and gcr_vals[0] > 0


def test_compute_poa_nonnegative():
    df = make_hourly_df()
    df = compute_solar_angles(df, lat_deg=10.0)
    df = compute_shading(df, beta_deg=15.0, d_row_m=4.0, H_m_m=2.0)
    out = compute_POA(df, beta_deg=15.0)
    assert "G_eff" in out.columns and "G_front" in out.columns
    assert out["G_eff"].min() >= 0
    assert out["G_front"].min() >= 0


def test_modules_per_ha_positive():
    n = modules_per_ha(beta_deg=20.0, d_row_m=4.0)
    assert isinstance(n, int)
    assert n > 0
