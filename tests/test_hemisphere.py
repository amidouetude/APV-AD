"""
Hemisphere orientation of the plane-of-array model.

Guards the 2026-09 correction in compute_POA: the collector faces the
equator, i.e. due south when lat >= 0 and due north when lat < 0. Before
the fix the due-south form was applied unconditionally, which cost up to
32 % of annual POA at 35 deg S (audit/phase0/).

The tests are written against CLOSED FORMS and an independent oracle
(pvlib, when installed) rather than stored expected numbers, so they stay
meaningful if parameters change and still fail if the physics breaks.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "notebooks"))

from simulation_functions import compute_POA, compute_solar_angles  # noqa: E402


# ----------------------------------------------------------------------
def _frame(dec_deg, ha_deg, dni=1000.0, ghi=0.0, dhi=0.0, f_shad=0.0):
    """Minimal frame carrying exactly what compute_POA reads."""
    n = len(np.atleast_1d(dec_deg))
    return pd.DataFrame({
        "declination_deg": np.atleast_1d(dec_deg).astype(float),
        "hour_angle_deg": np.atleast_1d(ha_deg).astype(float),
        "GHI": np.full(n, ghi, dtype=float),
        "DNI": np.full(n, dni, dtype=float),
        "DHI": np.full(n, dhi, dtype=float),
        "F_shad": np.full(n, f_shad, dtype=float),
    })


def _cos_aoi(lat, beta, dec, ha):
    """cos(theta_i) recovered from the beam term, with DNI = 1000."""
    out = compute_POA(_frame(dec, ha), beta, lat)
    return float(out["G_front"].iloc[0]) / 1000.0


# ======================================================================
# Closed-form identity at solar noon
# ======================================================================
@pytest.mark.parametrize("lat,beta,dec", [
    (37.87, 20.0, 0.0),      # Konya, equinox
    (12.37, 15.0, 10.0),     # Ouagadougou
    (47.99, 25.0, -20.0),    # Freiburg, winter
    (-15.42, 15.0, 0.0),     # Lusaka
    (-26.20, 20.0, 12.0),    # Johannesburg
    (-33.92, 30.0, -23.0),   # Cape Town, southern summer
])
def test_noon_matches_closed_form(lat, beta, dec):
    """
    At the solar-noon hour angle the Duffie & Beckman expression collapses
    to cos(phi_eff - delta), with phi_eff = phi - beta north / phi + beta
    south. This is the identity the pre-2026-09 code violated below the
    equator.
    """
    phi_eff = lat - (1.0 if lat >= 0 else -1.0) * beta
    expected = max(np.cos(np.radians(phi_eff - dec)), 0.0)
    assert _cos_aoi(lat, beta, dec, 0.0) == pytest.approx(expected, abs=1e-12)


# ======================================================================
# The collector must face the equator, not always south
# ======================================================================
@pytest.mark.parametrize("lat", [-5.0, -15.0, -25.0, -35.0])
def test_southern_tilt_beats_flat_at_local_noon(lat):
    """
    In the southern hemisphere, tilting an equator-facing collector must
    IMPROVE noon incidence over a horizontal panel at the equinox. The old
    due-south form made it worse, which is the whole bug.
    """
    flat = _cos_aoi(lat, 0.0, 0.0, 0.0)
    tilted = _cos_aoi(lat, abs(lat), 0.0, 0.0)
    assert tilted > flat
    assert tilted == pytest.approx(1.0, abs=1e-12)   # beta = |lat| faces the sun exactly


def test_hemispheres_are_mirror_symmetric():
    """
    A site at +phi at declination +delta must see exactly what a site at
    -phi sees at -delta, for the same tilt. Pure geometric symmetry: any
    asymmetry means an orientation assumption is still baked in.
    """
    for lat, beta, dec in [(20.0, 15.0, 10.0), (35.0, 25.0, -18.0), (8.0, 5.0, 3.0)]:
        north = _cos_aoi(lat, beta, dec, 30.0)
        south = _cos_aoi(-lat, beta, -dec, 30.0)
        assert north == pytest.approx(south, abs=1e-12)


# ======================================================================
# The equator: np.sign(0) == 0 must not silently flatten the array
# ======================================================================
def test_equator_keeps_the_tilt():
    """
    REGRESSION. The first proposed patch used np.sign(lat), and np.sign(0.0)
    is 0 -- which cancels beta entirely and treats the array as horizontal
    at exactly lat = 0. Caught by the phase-0 audit before it shipped.
    """
    beta = 20.0
    at_equator = _cos_aoi(0.0, beta, 0.0, 0.0)
    flat = _cos_aoi(0.0, 0.0, 0.0, 0.0)

    assert at_equator != pytest.approx(flat, abs=1e-9), (
        "tilt was dropped at lat = 0 -- np.sign(0) bug is back"
    )
    assert at_equator == pytest.approx(np.cos(np.radians(beta)), abs=1e-12)


def test_equator_convention_is_continuous_from_the_north():
    """
    lat >= 0 is documented as northern. The value at exactly 0 must
    therefore equal the limit approached from +epsilon, not from -epsilon.
    """
    beta = 20.0
    assert _cos_aoi(0.0, beta, 0.0, 0.0) == pytest.approx(
        _cos_aoi(1e-9, beta, 0.0, 0.0), abs=1e-9)


# ======================================================================
# Northern hemisphere must not move at all
# ======================================================================
@pytest.mark.parametrize("lat", [2.5, 12.37, 20.0, 37.87, 47.99])
def test_northern_unchanged_vs_legacy_form(lat):
    """
    Non-regression: north of the equator the corrected expression must be
    bit-identical to the legacy (phi - beta) form. The four published sites
    are all northern, so any drift here would invalidate the manuscript.
    """
    beta, dec, ha = 22.0, 12.0, 45.0
    legacy = np.clip(
        np.sin(np.radians(dec)) * np.sin(np.radians(lat - beta))
        + np.cos(np.radians(dec)) * np.cos(np.radians(ha)) * np.cos(np.radians(lat - beta)),
        0, 1)
    assert _cos_aoi(lat, beta, dec, ha) == pytest.approx(float(legacy), abs=1e-15)


# ======================================================================
# Independent oracle
# ======================================================================
def test_against_pvlib():
    """
    Cross-check against pvlib's own angle-of-incidence, fed with the same
    analytic solar position so the comparison isolates the surface
    geometry from the solar-position model.

    Convention trap handled explicitly: Duffie & Beckman measures the
    surface azimuth from SOUTH, pvlib from NORTH clockwise.
    """
    pvlib = pytest.importorskip("pvlib")

    rng = np.random.default_rng(7)
    for _ in range(300):
        lat = float(rng.uniform(-38, 38))
        beta = float(rng.uniform(0, 40))
        dec = float(rng.uniform(-23.45, 23.45))
        ha = float(rng.uniform(-75, 75))

        phi, d, w = map(np.radians, (lat, dec, ha))
        cos_z = np.sin(phi)*np.sin(d) + np.cos(phi)*np.cos(d)*np.cos(w)
        theta_z = np.degrees(np.arccos(np.clip(cos_z, -1, 1)))
        if theta_z > 80:
            continue
        sin_z = np.sin(np.radians(theta_z))
        cos_gs = np.clip((cos_z*np.sin(phi) - np.sin(d)) / (sin_z*np.cos(phi)), -1, 1)
        az_pvlib = np.degrees(np.sign(ha) * abs(np.arccos(cos_gs))) + 180.0

        surface_az = 180.0 if lat >= 0 else 0.0     # equator-facing
        ref = np.clip(np.cos(np.radians(pvlib.irradiance.aoi(
            surface_tilt=beta, surface_azimuth=surface_az,
            solar_zenith=theta_z, solar_azimuth=az_pvlib))), 0, 1)

        assert _cos_aoi(lat, beta, dec, ha) == pytest.approx(float(ref), abs=1e-9)


# ======================================================================
# End to end, through the real solar-angle model
# ======================================================================
def test_southern_annual_poa_exceeds_legacy():
    """
    Full-year sanity check through compute_solar_angles: at 25 deg S an
    equator-facing array must collect more than the legacy due-south form
    would have. Uses a synthetic clear-sky year -- no climate file needed.
    """
    idx = pd.date_range("2001-01-01", periods=8760, freq="h", tz=None)
    df = pd.DataFrame({"GHI": 0.0, "DNI": 1000.0, "DHI": 0.0}, index=idx)
    df = compute_solar_angles(df, -25.0, 28.0)
    df["F_shad"] = 0.0

    lat, beta = -25.0, 25.0
    corrected = compute_POA(df, beta, lat)["G_front"].sum()

    dec_r = np.radians(df["declination_deg"].values)
    ha_r = np.radians(df["hour_angle_deg"].values)
    legacy_cos = np.clip(
        np.sin(dec_r)*np.sin(np.radians(lat - beta))
        + np.cos(dec_r)*np.cos(ha_r)*np.cos(np.radians(lat - beta)), 0, 1)
    legacy = (1000.0 * legacy_cos).sum()

    assert corrected > legacy * 1.2, (
        f"corrected {corrected:.0f} vs legacy {legacy:.0f} -- expected a large gain"
    )
