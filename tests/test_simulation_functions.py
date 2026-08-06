# tests/test_simulation_functions.py
# =============================================================================
# Expanded pytest coverage for simulation_functions.py.
#
# Added in response to peer review (item F): the original test suite covered
# geometry/POA and eLER end-to-end via notebook unit tests, but had no
# pytest coverage for compute_PV, compute_crop_yield, or compute_AD in
# isolation, and no CI-friendly fixture data (tests required the full
# 8760-hour PVGIS TMY CSVs, which are too large/slow for a CI smoke test).
#
# Updated in response to a second review round: tests were previously
# centered on a single representative site (Konya). Site-dependent tests
# are now parametrized across all four sites (Konya, Almeria, Ouagadougou,
# Freiburg) via the `site` fixture, so each runs 4x with that site's real
# (lat, lon) and SITES[site] economic/agronomic parameters. This exercises
# the longitude-correction fix (fix #5) across its full range — Konya
# (lon 32.49E, the largest longitude and the site most affected by that
# bug) and Almeria/Ouagadougou (lon near 0E, the sites least affected) are
# both covered, rather than only the single site used previously.
#
# These tests use small SYNTHETIC fixtures (one week of hourly data,
# hand-constructed) instead of real PVGIS data, so they run in well under a
# second and don't require any external files. They check known physical
# relationships and bounds, not exact reference numbers, since the synthetic
# inputs are not real climate data.
#
# Run with: pytest tests/test_simulation_functions.py -v
# =============================================================================

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "notebooks"))

from config_00 import SITES, PARAMS, SCENARIOS
from simulation_functions import (
    compute_solar_angles,
    compute_shading,
    compute_POA,
    compute_PV,
    compute_crop_yield,
    compute_ET_reduction,
    compute_AD,
    modules_per_ha,
    reference_pv_kWh,
    eLER_weighted,
    eLER_objective,
)

ALL_SITES = ["Konya", "Almeria", "Ouagadougou", "Freiburg"]


# =============================================================================
# Fixtures — small synthetic hourly datasets (no external files needed)
# =============================================================================

@pytest.fixture(params=ALL_SITES)
def site(request):
    """Parametrizes dependent tests/fixtures across all four study sites."""
    return request.param


@pytest.fixture
def lat(site):
    return SITES[site]["lat"]


@pytest.fixture
def lon(site):
    return SITES[site]["lon"]


@pytest.fixture
def one_week_df():
    """
    One synthetic week (June solstice week) of hourly climate data.
    Values are illustrative, not real PVGIS data — tests check physical
    relationships and bounds, not exact reference numbers. Intentionally
    NOT site-parametrized itself (the synthetic climate profile is
    generic); it is paired with each site's real (lat, lon) and
    SITES[site] parameters by the `geometry_df` / `site` fixtures instead,
    so the same synthetic weather is exercised against every site's
    geometry and economics.
    """
    idx = pd.date_range("2024-06-18", periods=24 * 7, freq="h")
    n = len(idx)
    hour = idx.hour.values

    # Simple bell-shaped daytime GHI/DNI/DHI profile, zero at night
    daylight = np.clip(np.sin(np.pi * (hour - 5) / 14), 0, None)
    GHI = 800 * daylight
    DNI = 700 * daylight
    DHI = 150 * daylight

    T2m = 15 + 10 * daylight + np.random.RandomState(0).normal(0, 0.5, n)
    WS  = np.full(n, 3.0)

    df = pd.DataFrame(
        {"GHI": GHI, "DNI": DNI, "DHI": DHI, "T2m": T2m, "WS": WS},
        index=idx,
    )
    # Entire synthetic week counts as "in season" for these unit tests
    df["in_season"] = True
    return df


@pytest.fixture
def geometry_df(one_week_df, lat, lon):
    """one_week_df with solar angles + shading + POA already computed,
    for the current parametrized site's (lat, lon)."""
    df = compute_solar_angles(one_week_df, lat, lon)
    df = compute_shading(df, beta_deg=20.0, d_row_m=3.03, H_m_m=1.5)
    df = compute_POA(df, beta_deg=20.0, lat_deg=lat)
    return df


# =============================================================================
# compute_solar_angles — parametrized across all four sites
# =============================================================================

def test_solar_elevation_bounds(one_week_df, lat, lon):
    df = compute_solar_angles(one_week_df, lat, lon)
    assert df["solar_elevation_deg"].max() <= 90.0
    assert df["solar_elevation_deg"].min() >= -90.0
    # June week: all four sites (12-48N) should have strongly positive
    # elevation around local solar noon
    midday = df.between_time("10:00", "14:00")
    assert (midday["solar_elevation_deg"] > 20).any()


def test_solar_elevation_zero_at_night(one_week_df, lat, lon):
    df = compute_solar_angles(one_week_df, lat, lon)
    night = df.between_time("00:00", "02:00")
    assert (night["solar_elevation_deg"] < 10).all()


def test_solar_angles_requires_lon_deg(one_week_df, lat):
    """Regression test for the longitude-correction fix: the old 2-arg
    signature (no lon_deg) silently assumed solar noon at 12:00 UTC
    everywhere. Calling with only lat_deg must now raise a TypeError.
    Parametrized across sites so this is checked at every latitude used
    in the study, not just one."""
    with pytest.raises(TypeError):
        compute_solar_angles(one_week_df, lat)  # missing required lon_deg


def test_solar_noon_shifts_with_longitude(one_week_df, site):
    """East of the Prime Meridian, true solar noon occurs BEFORE 12:00 UTC.
    Checked using each site's own real longitude against a lon=0 control
    at the same latitude, so the sign/magnitude of the shift is verified
    for every site in the study rather than only the most-affected one."""
    real_lon = SITES[site]["lon"]
    lat = SITES[site]["lat"]
    df_real = compute_solar_angles(one_week_df, lat, lon_deg=real_lon)
    df_greenwich = compute_solar_angles(one_week_df, lat, lon_deg=0.0)

    def peak_hour(df):
        day1 = df.loc[df.index.date == df.index.date[0]]
        return day1["solar_elevation_deg"].idxmax().hour

    if abs(real_lon) < 0.5:
        pytest.skip(f"{site}'s longitude ({real_lon}) is too close to 0 "
                    "for this test to meaningfully distinguish a shift")
    elif real_lon > 0:
        assert peak_hour(df_real) <= peak_hour(df_greenwich)
    else:
        assert peak_hour(df_real) >= peak_hour(df_greenwich)


# =============================================================================
# compute_POA — regression coverage for the azimuth fix, all four sites
# =============================================================================

def test_poa_requires_lat_deg(one_week_df, lat, lon):
    """compute_POA's signature must include lat_deg (regression test for the
    angle-of-incidence fix — calling with the old 2-arg signature should
    raise a TypeError, not silently ignore latitude)."""
    df = compute_solar_angles(one_week_df, lat, lon)
    df = compute_shading(df, 20.0, 3.03, 1.5)
    with pytest.raises(TypeError):
        compute_POA(df, 20.0)  # missing required lat_deg


def test_poa_nonnegative(geometry_df):
    assert (geometry_df["G_front"] >= 0).all()
    assert (geometry_df["G_rear"] >= 0).all()
    assert (geometry_df["G_eff"] >= 0).all()


def test_poa_zero_at_night(geometry_df):
    night = geometry_df.between_time("00:00", "02:00")
    # front POA should be ~0 when GHI/DNI are 0
    assert night["G_front"].max() < 1.0


# =============================================================================
# compute_PV — parametrized across all four sites
# =============================================================================

def test_compute_pv_nonnegative_energy(geometry_df):
    n_mod = modules_per_ha(20.0, 3.03)
    df = compute_PV(geometry_df, n_mod, f_PV_heat=0.05)
    assert (df["E_PV_kWh"] >= 0).all()
    assert (df["E_PV_heat_kWh"] >= 0).all()
    assert (df["E_PV_sold_kWh"] >= 0).all()


def test_compute_pv_heat_sold_sum_to_total(geometry_df):
    n_mod = modules_per_ha(20.0, 3.03)
    df = compute_PV(geometry_df, n_mod, f_PV_heat=0.20)
    total = df["E_PV_heat_kWh"] + df["E_PV_sold_kWh"]
    np.testing.assert_allclose(total.values, df["E_PV_kWh"].values, rtol=1e-9)


def test_compute_pv_zero_at_night(geometry_df):
    n_mod = modules_per_ha(20.0, 3.03)
    df = compute_PV(geometry_df, n_mod, f_PV_heat=0.05)
    night = df.between_time("00:00", "02:00")
    assert night["E_PV_kWh"].max() < 1.0


def test_compute_pv_scales_with_module_count(geometry_df):
    df_a = compute_PV(geometry_df, n_modules=100, f_PV_heat=0.05)
    df_b = compute_PV(geometry_df, n_modules=200, f_PV_heat=0.05)
    # double the modules -> ~double the energy (temperature term is per-module,
    # not per-array, so this should hold closely)
    ratio = df_b["E_PV_kWh"].sum() / max(df_a["E_PV_kWh"].sum(), 1e-9)
    assert 1.9 < ratio < 2.1


# =============================================================================
# compute_ET_reduction — parametrized across all four sites
# =============================================================================

def test_et_reduction_physically_plausible_daily_rate(one_week_df, site, lat, lon):
    """Regression test for the missing-0.408-factor bug: ET0 must fall in a
    physically plausible daily range (roughly 1-10 mm/day for any climate;
    the bug produced ~11+ mm/day at a temperate semi-arid site, which is
    implausible). Parametrized across all four sites, spanning the full
    12-48N latitude range used in the study."""
    df = compute_solar_angles(one_week_df, lat, lon)
    df = compute_shading(df, beta_deg=20.0, d_row_m=3.03, H_m_m=1.5)
    result = compute_ET_reduction(df, site)
    days = one_week_df["in_season"].sum() / 24
    daily_rate = result["ET0_season_mm"] / days
    assert 0.5 < daily_rate < 10.0, (
        f"[{site}] ET0 daily rate {daily_rate:.2f} mm/day is not physically "
        "plausible — check the Hargreaves-Samani MJ-to-mm conversion factor"
    )


# =============================================================================
# compute_crop_yield — parametrized across all four sites
# =============================================================================

def test_crop_yield_bounded(geometry_df):
    result = compute_crop_yield(geometry_df, PAR_sat=174.0)
    assert 0.0 <= result["LER_crop"] <= 1.1  # allow tiny numerical slack


def test_crop_yield_no_shading_is_higher_than_with_shading(one_week_df, lat, lon):
    df = compute_solar_angles(one_week_df, lat, lon)
    df_noshade = compute_shading(df.copy(), beta_deg=20.0, d_row_m=20.0, H_m_m=0.01)
    df_shade   = compute_shading(df.copy(), beta_deg=20.0, d_row_m=3.03,  H_m_m=1.5)

    r_noshade = compute_crop_yield(df_noshade, PAR_sat=174.0)
    r_shade   = compute_crop_yield(df_shade,   PAR_sat=174.0)

    # More shading should never increase LER_crop
    assert r_shade["LER_crop"] <= r_noshade["LER_crop"] + 1e-9


def test_crop_yield_zero_ghi_returns_zero():
    # Not site-dependent (pure zero-input edge case) — left unparametrized.
    idx = pd.date_range("2024-06-18", periods=24, freq="h")
    df = pd.DataFrame({
        "GHI": np.zeros(24), "F_shad": np.zeros(24), "in_season": True,
    }, index=idx)
    result = compute_crop_yield(df, PAR_sat=174.0)
    assert result["LER_crop"] == 0.0


# =============================================================================
# compute_AD — parametrized across all four sites, including AD_overrides
# =============================================================================

def test_compute_ad_returns_expected_keys(one_week_df, site):
    result = compute_AD(one_week_df, site, V_dig_m3=4.0, HRT_days=25.0,
                        LER_crop=0.37, scenario="S0")
    for key in ["BMP_avg_NmL_gVS", "VS_input_kg_ha", "OLR_kgVS_m3d",
                "E_biogas_kWh_ha", "Q_heat_kWh_ha"]:
        assert key in result


def test_compute_ad_nonnegative(one_week_df, site):
    result = compute_AD(one_week_df, site, V_dig_m3=4.0, HRT_days=25.0,
                        LER_crop=0.37, scenario="S0")
    assert result["E_biogas_kWh_ha"] >= 0
    assert result["VS_input_kg_ha"] >= 0
    assert result["BMP_avg_NmL_gVS"] >= 0


def test_compute_ad_overrides_default_matches_no_override(one_week_df, site):
    """AD_overrides=None must reproduce the exact same result as passing no
    argument at all (backward-compatibility regression test)."""
    r1 = compute_AD(one_week_df, site, 4.0, 25.0, 0.37, "S0")
    r2 = compute_AD(one_week_df, site, 4.0, 25.0, 0.37, "S0", AD_overrides=None)
    assert r1 == r2


def test_compute_ad_overrides_changes_output(one_week_df, site):
    baseline = compute_AD(one_week_df, site, 4.0, 25.0, 0.37, "S0")
    perturbed = compute_AD(one_week_df, site, 4.0, 25.0, 0.37, "S0",
                           AD_overrides={"BMP_ref": 200.0})
    assert perturbed["E_biogas_kWh_ha"] < baseline["E_biogas_kWh_ha"]


def test_compute_ad_higher_capture_efficiency_increases_biogas(one_week_df, site):
    low = compute_AD(one_week_df, site, 4.0, 25.0, 0.37, "S0",
                     AD_overrides={"eta_cap": 0.70})
    high = compute_AD(one_week_df, site, 4.0, 25.0, 0.37, "S0",
                      AD_overrides={"eta_cap": 1.00})
    assert high["E_biogas_kWh_ha"] > low["E_biogas_kWh_ha"]


# =============================================================================
# eLER_weighted — pure function of its numeric arguments, not site-dependent
# =============================================================================

def test_eler_weighted_equal_weights_reproduces_sum():
    val = eLER_weighted(0.4, 0.6, 0.5, w_crop=1.0, w_PV=1.0, w_bio=1.0)
    assert val == pytest.approx(0.4 + 0.6 + 0.5, rel=1e-9)


def test_eler_weighted_normalizes_arbitrary_scale():
    # Same relative weights at two different absolute scales must give the
    # same result, since weights are renormalized to sum to 3 internally.
    v1 = eLER_weighted(0.4, 0.6, 0.5, w_crop=1.0, w_PV=2.0, w_bio=1.0)
    v2 = eLER_weighted(0.4, 0.6, 0.5, w_crop=10.0, w_PV=20.0, w_bio=10.0)
    assert v1 == pytest.approx(v2, rel=1e-9)


def test_eler_weighted_rejects_nonpositive_weight_sum():
    with pytest.raises(ValueError):
        eLER_weighted(0.4, 0.6, 0.5, w_crop=-1.0, w_PV=-1.0, w_bio=-1.0)


# =============================================================================
# eLER_objective — end-to-end smoke test on synthetic data, all four sites
# =============================================================================

def test_eler_objective_runs_end_to_end(one_week_df, site):
    x = np.array([25.0, 3.03, 1.5, 3.5, 25.0, 0.05])
    # NOTE: beta=25 and V_dig=3.5 (not round numbers like 20/4.0) are
    # deliberately chosen to be feasible at every one of the four
    # sites' differing min_beta_deg and OLR constraints (Freiburg's
    # lower fruit_yield_t_ha, e.g., makes OLR infeasible at V_dig=4.0).
    # This was caught by site-parametrization itself: a design point
    # that worked for 3/4 sites by coincidence was not actually
    # robust across the full parameter space this study covers.
    result = eLER_objective(x, one_week_df, site, return_full=True)
    assert result is not None
    assert result["eLER_defA"] > 0
    assert result["eLER_weighted"] == pytest.approx(result["eLER_defA"], rel=1e-9)


def test_eler_objective_infeasible_tilt_returns_penalty(one_week_df, site):
    x = np.array([5.0, 3.03, 1.5, 4.0, 25.0, 0.05])  # below min_beta_deg
    val = eLER_objective(x, one_week_df, site)
    assert val < 0  # infeasibility penalty


def test_eler_objective_weighted_matches_manual_calc(one_week_df, site):
    x = np.array([25.0, 3.03, 1.5, 3.5, 25.0, 0.05])  # see note above
    w = (0.5, 2.0, 0.5)
    result = eLER_objective(x, one_week_df, site, return_full=True,
                            w_crop=w[0], w_PV=w[1], w_bio=w[2])
    # Recompute from the *rounded* component values returned in the dict —
    # allow a small absolute tolerance since eLER_objective rounds each
    # component to 3 decimals before this reconstruction, while its own
    # internal eLER_weighted value is computed from unrounded components.
    expected = eLER_weighted(result["LER_crop"], result["LER_PV_defA"],
                             result["LER_biogas"], *w)
    assert result["eLER_weighted"] == pytest.approx(expected, abs=2e-3)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
