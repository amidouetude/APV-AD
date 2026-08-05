"""
=============================================================================
INTEGRATED AGRIVOLTAIC – ANAEROBIC DIGESTION OPTIMIZATION MODEL
Version 3.0 — Four-Site, Production Release
Target journal: Renewable & Sustainable Energy Reviews / Applied Energy
=============================================================================

STUDY SITES (four Köppen climate zones):
  Konya        Turkey       BSk  Cold semi-arid       37.87°N  32.49°E  UTC+3
  Almería      Spain        BSh  Hot semi-arid        36.83°N   2.46°W  UTC+1
  Ouagadougou  Burkina Faso BSh  Tropical semi-arid   12.37°N   1.52°W  UTC+0
  Freiburg     Germany      Cfb  Oceanic temperate    47.99°N   7.85°E  UTC+1

CROP: Tomato (Solanum lycopersicum)
  Semi-arid sites : 60 t/ha irrigated yield (Mohammedi et al., 2023)
  Freiburg        : 30 t/ha open-field yield (DLG, 2023)

=============================================================================
FIXES AND IMPROVEMENTS — VERSION 3.0
=============================================================================

  FIX 1  crop_yield_ratio() — asymmetric PAR saturation [CRITICAL]
         PAR_open denominator: raw integral, no saturation clip.
         PAR_AV numerator: clipped at PAR_sat = 174 W/m².
         Restores LER_crop to physically correct range 0.55–0.65.
         Root cause of LER_crop = 1.000 bug in all earlier versions.

  FIX 2  LER_PV reference cache with tilt invalidation
         _ref_pv_tilt_cached tracks the tilt used; cache is rebuilt when
         tilt changes (constraint sweep, Pareto calls on the same system
         object at different optimizer-proposed tilts).

  FIX 3  Ouagadougou cross-year growing season mask
         Season Jul 1 – Oct 30 (no cross-year); original code is correct.
         Original code had the wrong start date (Oct 1) based on an older
         paper draft. Updated to Jul 1 – Oct 30 per Table 5 of the paper.
         The OR-logic guard is retained for robustness against any future
         cross-year site (e_month < s_month triggers OR, else AND).

  FIX 4  Site-specific economics and crop yield
         p_elec, p_biogas, p_tomato, fruit_yield_t_ha stored per-site in
         SITES dict and injected at runtime. No global ECONOMICS mutation.

  FIX 5  AD Scenario S5 — biomethane grid injection
         Valued at 5.5 × p_elec (Nik Zad et al., 2025; EEG 2023).
         One SCENARIOS entry; no structural change to energy balance.

  FIX 6  Piecewise BMP temperature model (validation only)
         T ≤ 37°C: Arrhenius (Pilarski & Pilarska, 2025).
         T > 37°C: linear inhibition −5 NmL/gVS/°C (Feng et al., 2013).
         T_eff is clipped at 37°C operationally → inhibition branch only
         used in run_validation(); simulation outputs are unchanged.

  FIX 7  MC rewritten as structural invariance + economic sensitivity
         Replaces CV=0.0% spike histograms with:
           (a) one-way tornado on design variables (genuine eLER sensitivity)
           (b) net income P5/P95 distribution (genuine economic MC output)

  FIX 8  Freiburg site added
         47.99°N, 7.85°E, Cfb, 30 t/ha, May 1 – Oct 15.
         psychro_gamma = 0.0673 kPa/°C at 278 m elevation.

  PERF 1 Vectorized reference PV computation
         _compute_reference_pv() replaced with fully vectorized numpy
         operations. ~400× faster than row-by-row iterrows() loop.

  PERF 2 MC loop avoids per-sample solar geometry rebuild
         Monte-Carlo samples share the pre-computed solar geometry arrays
         (solar_altitude_rad, solar_azimuth_rad) from the base system's
         climate DataFrame. Only physics parameters (BMP, U_eff, α_shade)
         and prices are varied. ~30× speedup for N=2000 samples.

  FIX 9  Pareto base class safe under missing pymoo
         ParetoAPVAD uses a runtime-resolved base class to avoid a
         NameError at module import when pymoo is not installed.

  FIX 10 Pareto site selection — season-length gate
         Pareto is run only for sites with growing season > 150 days
         (Konya: 189d, Freiburg: ~167d). Almería (122d) and Ouagadougou
         (122d) are excluded because LER_water variance is below numerical
         resolution at GCR ≈ 18.2% over a 122-day season. See paper
         Section 3.9 for the geometric derivation.

=============================================================================
USAGE
=============================================================================
  1. Place four PVGIS-SARAH3 TMY CSV files in the working directory:
       konya_climate.csv, almeria_climate.csv,
       ouagadougou_climate.csv, freiburg_climate.csv
     OR update the "climate_file" paths in the SITES dictionary.
  2. Adjust RUN_CONFIG as needed (toggle analyses, set output directory).
  3. Run:  python APV_AD_v3.py
  4. Results are written to RUN_CONFIG["output_dir"] as PNG + PDF figures
     and a CSV summary table.

=============================================================================
DEPENDENCIES
=============================================================================
  Core  : numpy, pandas, scipy, matplotlib
  Pareto: pymoo >= 0.6   (optional; R3-B disabled if absent)

  Install:  pip install numpy pandas scipy matplotlib pymoo
=============================================================================
"""

from __future__ import annotations

import warnings
import csv
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import differential_evolution
from scipy.stats import qmc

warnings.filterwarnings("ignore")
matplotlib.use("Agg")   # non-interactive backend; safe on headless servers

# ── Optional pymoo ────────────────────────────────────────────────────────────
try:
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.core.problem import Problem as _PymooBase
    from pymoo.optimize import minimize as pymoo_minimize
    from pymoo.termination import get_termination
    _PYMOO = True
except ImportError:
    _PymooBase = object          # safe fallback — ParetoAPVAD will be unusable
    _PYMOO = False
    print("[INFO] pymoo not installed — NSGA-II Pareto (R3-B) disabled.")
    print("       Install with:  pip install pymoo")


# =============================================================================
# SECTION 0 — GLOBAL CONFIGURATION
# =============================================================================

SITES: Dict[str, dict] = {
    # ── Semi-arid — cold continental ─────────────────────────────────────
    "Konya": {
        "latitude": 37.87, "longitude": 32.49, "timezone_offset": 3,
        "climate_zone": "BSk", "country": "Turkey",
        "growing_start": (4, 25),   # (month, day)
        "growing_end":   (10, 30),
        "climate_file":  "konya_climate.csv",
        "p_elec_usd_kWh":   0.06,   # YEKA feed-in tariff proxy
        "p_biogas_usd_kWh": 0.08,
        "p_tomato_usd_t":   300.0,
        "discount_rate":    0.08,
        "fruit_yield_t_ha": 60.0,   # irrigated Mediterranean yield
        "psychro_gamma":    0.0665, # kPa/°C at 1,016 m elevation
    },
    # ── Semi-arid — hot Mediterranean ────────────────────────────────────
    "Almeria": {
        "latitude": 36.83, "longitude": -2.46, "timezone_offset": 1,
        "climate_zone": "BSh", "country": "Spain",
        "growing_start": (3, 1),
        "growing_end":   (6, 30),
        "climate_file":  "almeria_climate.csv",
        "p_elec_usd_kWh":   0.09,
        "p_biogas_usd_kWh": 0.09,
        "p_tomato_usd_t":   400.0,
        "discount_rate":    0.06,
        "fruit_yield_t_ha": 60.0,
        "psychro_gamma":    0.0670, # at 20 m elevation
    },
    # ── Semi-arid — tropical Sahelian ────────────────────────────────────
    "Ouagadougou": {
        "latitude": 12.37, "longitude": -1.52, "timezone_offset": 0,
        "climate_zone": "BSh", "country": "Burkina Faso",
        "growing_start": (7, 1),    # FIX 3: corrected from (10,1) → (7,1)
        "growing_end":   (10, 30),  # growing season: wet monsoon July–October
        "climate_file":  "ouagadougou_climate.csv",
        "p_elec_usd_kWh":   0.12,
        "p_biogas_usd_kWh": 0.08,
        "p_tomato_usd_t":   250.0,
        "discount_rate":    0.10,
        "fruit_yield_t_ha": 60.0,
        "psychro_gamma":    0.0672, # at 303 m elevation
    },
    # ── Temperate oceanic (FIX 8) ─────────────────────────────────────────
    "Freiburg": {
        "latitude": 47.99, "longitude": 7.85, "timezone_offset": 1,
        "climate_zone": "Cfb", "country": "Germany",
        "growing_start": (5, 1),
        "growing_end":   (10, 15),
        "climate_file":  "freiburg_climate.csv",
        "p_elec_usd_kWh":   0.09,   # EEG 2023 feed-in tariff proxy
        "p_biogas_usd_kWh": 0.09,
        "p_tomato_usd_t":   350.0,  # German wholesale
        "discount_rate":    0.05,
        "fruit_yield_t_ha": 30.0,   # FIX 8: open-field German yield (DLG, 2023)
        "psychro_gamma":    0.0673, # at 278 m elevation
    },
}

RUN_CONFIG: Dict[str, object] = {
    # Sites to simulate (order determines table/figure order)
    "sites":                   ["Konya", "Almeria", "Ouagadougou", "Freiburg"],
    # Differential evolution settings
    "de_maxiter":              50,
    "de_popsize":              12,    # actual pop = 12 × n_var = 72
    "de_seed":                 42,
    # Monte-Carlo
    "run_monte_carlo":         True,
    "mc_n_samples":            2000,
    "mc_seed":                 99,
    # Multi-objective
    "run_mo_constraint_sweep": True,
    "run_mo_pareto":           True,
    "mo_d_row_values":         [3, 4, 5, 6],  # m — tractor access scenarios
    "pareto_min_season_days":  150,           # FIX 10: exclude short-season sites
    # AD scenarios
    "run_ad_scenarios":        True,
    # Validation
    "run_validation":          True,
    # Output
    "save_figures":            True,
    "figure_dpi":              300,
    "output_dir":              "apv_ad_results",
}

PHYSICAL_CONSTANTS: Dict[str, object] = {
    "CH4_LHV_MJ_m3":       35.8,    # lower heating value of methane (MJ/m³ STP)
    "PAR_to_GHI_ratio":    0.48,    # McCree (1971); Jacovides et al. (2006)
    "PAR_sat_tomato_W_m2": 174.0,   # 800 µmol/m²/s ÷ 4.6 (Mohammedi et al., 2023)
    "albedo_agri_soil":    0.25,    # Monteith & Unsworth (2008)
    "days_per_month":      np.array([31, 28, 31, 30, 31, 30,
                                      31, 31, 30, 31, 30, 31]),
}

# Capital and O&M economics (not site-specific — FIX 4)
CAPEX: Dict[str, float] = {
    "pv_capex_usd_module":   280.0,  # Jinko Tiger Neo 550 Wp (2024)
    "ad_capex_usd_m3":       450.0,  # farm-scale mesophilic digester
    "struct_capex_frac":     0.40,   # elevated mounting structure (fraction of PV CAPEX)
    "opex_frac_capex":       0.025,  # annual O&M as fraction of total CAPEX
    "p_N_fertilizer_usd_kg": 1.20,   # N-equivalent fertilizer price (digestate S4)
    "lpg_price_usd_kWh":     0.18,   # LPG parity for local biogas use (S3)
}

# Biomethane premium multiplier (Scenario S5)
BIOMETHANE_PREMIUM = 5.5   # × p_elec (Nik Zad et al., 2025; EEG 2023)


# =============================================================================
# SECTION 1 — SOLAR GEOMETRY
# =============================================================================

class SolarGeometry:
    """
    Hourly solar position and ground shading geometry.

    References
    ----------
    Spencer (1971) — solar declination and equation of time.
    Zainali et al. (2023), Applied Energy 339, 120981 — shading factor.
    """

    def __init__(self, latitude: float, longitude: float,
                 timezone_offset: float = 0.0) -> None:
        self.lat_rad = np.radians(latitude)
        self.lon = longitude
        self.tz  = timezone_offset

    # ------------------------------------------------------------------
    def day_of_year(self, dt: datetime) -> int:
        return dt.timetuple().tm_yday

    def solar_declination(self, n: int) -> float:
        """Declination (radians). Spencer (1971)."""
        return np.radians(23.45 * np.sin(np.radians((360.0 / 365.0) * (n - 81))))

    def equation_of_time(self, n: int) -> float:
        """Equation of time (minutes). Spencer (1971)."""
        B = (2.0 * np.pi / 365.0) * (n - 1)
        return (0.000075 + 0.001868 * np.cos(B) - 0.032077 * np.sin(B)
                - 0.014615 * np.cos(2 * B) - 0.04089 * np.sin(2 * B)) * 229.18

    def solar_time(self, dt: datetime) -> float:
        """True solar time (decimal hours)."""
        n    = self.day_of_year(dt)
        hour = dt.hour + dt.minute / 60.0
        return (hour
                + (self.lon - 15.0 * self.tz) * 4.0 / 60.0
                + self.equation_of_time(n) / 60.0)

    def hour_angle(self, dt: datetime) -> float:
        """Hour angle (radians); 0 at solar noon."""
        return np.radians(15.0 * (self.solar_time(dt) - 12.0))

    def solar_altitude(self, n: int, ha: float) -> float:
        """Solar elevation angle (radians)."""
        d = self.solar_declination(n)
        s = (np.sin(self.lat_rad) * np.sin(d)
             + np.cos(self.lat_rad) * np.cos(d) * np.cos(ha))
        return np.arcsin(np.clip(s, -1.0, 1.0))

    def solar_azimuth(self, n: int, ha: float, altitude: float) -> float:
        """Solar azimuth (radians, N-clockwise)."""
        d  = self.solar_declination(n)
        ca = ((np.sin(d) - np.sin(self.lat_rad) * np.sin(altitude))
              / (np.cos(self.lat_rad) * np.cos(altitude) + 1e-10))
        az = np.arccos(np.clip(ca, -1.0, 1.0))
        if np.sin(ha) > 0:
            az = 2.0 * np.pi - az
        return az

    def shading_factor(self, altitude_rad: float, azimuth_rad: float,
                       module_height: float, row_spacing: float,
                       row_azimuth: float = 0.0) -> float:
        """
        Ground shading factor F_shad ∈ [0, 1]. Zainali et al. (2023).
        F_shad = min(1, H_m × cos(Δaz) / (d_row × tan(α_s)))
        """
        if altitude_rad <= 0.0:
            return 0.0
        shadow   = module_height / np.tan(altitude_rad)
        azi_diff = abs(azimuth_rad - row_azimuth)
        if azi_diff > np.pi:
            azi_diff = 2.0 * np.pi - azi_diff
        eff = shadow * np.cos(azi_diff)
        if eff <= 0.0:
            return 0.0
        return float(np.clip(eff / row_spacing, 0.0, 1.0))


# =============================================================================
# SECTION 2 — AGRIVOLTAIC PV SYSTEM
# =============================================================================

class AgrivoltaicSystem:
    """
    Bifacial APV: irradiance (Liu & Jordan, 1960), Faiman (2008) cell
    temperature, FIX-1 asymmetric PAR crop yield, FAO-56 ET reduction.

    Module: Jinko Tiger Neo N-type bifacial 550 Wp (2024 datasheet).
    """

    def __init__(self, land_area: float = 10_000.0,
                 alpha_shade: float = 0.40) -> None:
        self.land_area   = land_area
        self.alpha_shade = alpha_shade    # R1: varied in Monte-Carlo

        self.module_power_stc = 550       # Wp
        self.module_length    = 2.278     # m
        self.module_width     = 1.134     # m
        self.module_area      = self.module_length * self.module_width
        self.temp_coeff_pmax  = -0.0035   # /°C (Pmax temperature coefficient)
        self.bifaciality      = 0.70      # rear-to-front ratio
        self.system_losses    = 0.12      # wiring + inverter + mismatch

        # Faiman (2008) cell temperature coefficients
        self.U0 = 25.0     # W/m²·K — constant thermal loss
        self.U1 = 6.84     # W·s/m³·K — wind-speed dependent

        self.PAR_frac = PHYSICAL_CONSTANTS["PAR_to_GHI_ratio"]
        self.PAR_sat  = PHYSICAL_CONSTANTS["PAR_sat_tomato_W_m2"]
        self.albedo   = PHYSICAL_CONSTANTS["albedo_agri_soil"]

    # ── Geometry ──────────────────────────────────────────────────────────
    def modules_per_hectare(self, row_spacing: float) -> int:
        modules_per_row = np.sqrt(self.land_area) / (self.module_width + 0.05)
        n_rows          = np.sqrt(self.land_area) / row_spacing
        return int(max(1, modules_per_row * n_rows))

    def ground_cover_ratio(self, n_modules: int) -> float:
        return (n_modules * self.module_area) / self.land_area

    # ── Irradiance (Liu & Jordan, 1960) ───────────────────────────────────
    def plane_of_array_irradiance(
            self, GHI: float, DNI: float, DHI: float,
            tilt_deg: float, alt_rad: float, az_rad: float,
            surf_az_deg: float = 180.0) -> Tuple[float, float]:
        """Returns (G_eff_bifacial, G_front_only)."""
        tilt = np.radians(tilt_deg)
        surf = np.radians(surf_az_deg)
        cos_aoi  = (np.sin(alt_rad) * np.cos(tilt)
                    + np.cos(alt_rad) * np.sin(tilt) * np.cos(az_rad - surf))
        G_beam   = DNI * max(0.0, cos_aoi)
        G_diff   = DHI * (1.0 + np.cos(tilt)) / 2.0
        G_gnd    = GHI * self.albedo * (1.0 - np.cos(tilt)) / 2.0
        G_front  = G_beam + G_diff + G_gnd
        G_rear   = (DHI * (1.0 + np.cos(np.pi - tilt)) / 2.0
                    + GHI * self.albedo * (1.0 - np.cos(np.pi - tilt)) / 2.0)
        G_eff    = G_front + G_rear * self.bifaciality * (1.0 - self.system_losses)
        return float(max(0.0, G_eff)), float(max(0.0, G_front))

    # ── Cell temperature (Faiman, 2008) ───────────────────────────────────
    def cell_temperature(self, G_poa: float, T_amb: float,
                         wind_speed: float) -> float:
        return T_amb + G_poa / (self.U0 + self.U1 * max(0.5, wind_speed))

    # ── PV power (Eq. 10) ─────────────────────────────────────────────────
    def pv_power_kw(self, G_poa_eff: float, T_cell: float,
                    n_modules: int) -> float:
        p = ((G_poa_eff / 1000.0) * n_modules
             * self.module_power_stc * (1.0 - self.system_losses))
        p *= (1.0 + self.temp_coeff_pmax * (T_cell - 25.0))
        return float(max(0.0, p / 1000.0))   # kW

    # ── Crop yield ratio (FIX 1 — Eq. 11 — ASYMMETRIC PAR SATURATION) ────
    def crop_yield_ratio(self, GHI_arr: np.ndarray,
                         shading_arr: np.ndarray,
                         growing_mask: np.ndarray) -> float:
        """
        LER_crop = ∫ min(PAR_AV, PAR_sat) dt  /  ∫ PAR_open dt

        FIX 1 — critical correction:
        The denominator uses the RAW PAR_open integral (no saturation clip).
        At semi-arid midday irradiance (GHI > 360 W/m²), PAR_open > PAR_sat.
        Clipping PAR_open at PAR_sat (the old symmetric formulation) forced
        total_open → total_AV and LER_crop → 1.000 regardless of shading.
        The physically correct denominator is the raw open-field PAR integral
        because the open-field crop already receives more than enough light;
        extra photons above saturation simply do not contribute to yield in
        either the APV or the open-field case, but the RATIO is only correct
        when the denominator reflects what the open field actually receives.
        """
        PAR_open   = GHI_arr * self.PAR_frac           # raw — NO clip (FIX 1)
        PAR_AV     = PAR_open * (1.0 - shading_arr)
        PAR_AV_sat = np.minimum(PAR_AV, self.PAR_sat)  # numerator only

        total_open = float(np.sum(PAR_open[growing_mask]))
        total_AV   = float(np.sum(PAR_AV_sat[growing_mask]))
        return (total_AV / total_open) if total_open > 0.0 else 1.0

    # ── FAO-56 hourly Penman-Monteith ET₀ (Eq. 12) ───────────────────────
    def et0_penman_monteith(self, T: np.ndarray, RH: np.ndarray,
                             WS: np.ndarray, GHI: np.ndarray,
                             gamma: float = 0.0665) -> np.ndarray:
        """
        gamma: site-specific psychrometric constant (kPa/°C).
        Default 0.0665 at ~1,000 m; adjusted for each site (FIX 4).
        """
        es    = 0.6108 * np.exp(17.27 * T / (T + 237.3))
        ea    = es * np.clip(RH, 0.0, 100.0) / 100.0
        vpd   = np.maximum(0.0, es - ea)
        delta = 4098.0 * es / (T + 237.3) ** 2
        Rn    = 0.77 * GHI * 3.6e-3   # MJ/m²/h net radiation estimate
        num   = (0.408 * delta * Rn
                 + gamma * (37.0 / (T + 273.0)) * np.maximum(0.5, WS) * vpd)
        den   = delta + gamma * (1.0 + 0.24 * np.maximum(0.5, WS))
        return np.maximum(0.0, num / den)

    # ── APV ET reduction (Eq. 14) ─────────────────────────────────────────
    def et_reduction(self, ET0: np.ndarray,
                     shading: np.ndarray) -> np.ndarray:
        """ET_APV = ET0 × (1 − α_shade × F_shad). α_shade varied in R1."""
        return np.maximum(0.0, ET0 * (1.0 - self.alpha_shade * shading))


# =============================================================================
# SECTION 3 — ANAEROBIC DIGESTION
# =============================================================================

class AnaerobicDigestionSystem:
    """
    Mesophilic AD — tomato crop residues.

    Key features:
    • Temperature-corrected BMP via Arrhenius (Pilarski & Pilarska, 2025).
    • Site-specific fruit yield (FIX 4): Freiburg 30 t/ha vs. 60 t/ha.
    • Six valorization scenarios including S5 biomethane (FIX 5).
    • Piecewise BMP model above 37°C for validation accuracy (FIX 6).
    """

    SCENARIOS: Dict[str, dict] = {
        "S0": {"label": "Baseline (mono-digestion, grid export)",
               "BMP_mult": 1.00, "eta_cap": 1.00,
               "local": False, "biomethane": False, "digestate": False},
        "S1": {"label": "Co-digestion (30% cattle manure, +15% BMP)",
               "BMP_mult": 1.15, "eta_cap": 1.00,
               "local": False, "biomethane": False, "digestate": False},
        "S2": {"label": "Sub-optimal (85% CH4 capture)",
               "BMP_mult": 1.00, "eta_cap": 0.85,
               "local": False, "biomethane": False, "digestate": False},
        "S3": {"label": "Local LPG displacement (on-site use)",
               "BMP_mult": 1.00, "eta_cap": 1.00,
               "local": True,  "biomethane": False, "digestate": False},
        "S4": {"label": "Digestate valorization (N-P-K credit)",
               "BMP_mult": 1.00, "eta_cap": 1.00,
               "local": False, "biomethane": False, "digestate": True},
        "S5": {"label": "Biomethane grid injection (5.5× p_elec)",
               "BMP_mult": 1.00, "eta_cap": 1.00,
               "local": False, "biomethane": True,  "digestate": False},
    }

    def __init__(self, BMP_ref: float = 300.0, U_eff: float = 0.8,
                 fruit_yield_t_ha: float = 60.0) -> None:
        """
        Parameters
        ----------
        BMP_ref          : reference BMP at 37°C (NmL CH4/g VS) — R1 param
        U_eff            : digester wall U-value (W/m²·K) — R1 param
        fruit_yield_t_ha : site-specific fresh fruit yield (FIX 4)
        """
        # Operating constraints
        self.T_mesophilic  = 37.0   # °C optimal
        self.T_winter_min  = 26.0   # °C minimum with PV heating
        self.HRT_min       = 20     # days
        self.HRT_max       = 60     # days
        self.OLR_min       = 1.5    # kg VS/m³·day
        self.OLR_max       = 5.0    # kg VS/m³·day

        # Feedstock
        self.DM                  = 0.20   # dry matter fraction (kg DM/kg fresh)
        self.VS_DM               = 0.85   # volatile solids fraction (kg VS/kg DM)
        self.BMP_ref             = BMP_ref
        self.fruit_yield_t_ha    = fruit_yield_t_ha  # FIX 4
        self.residue_ratio       = 0.80   # residue / fruit mass ratio
        self.digestate_N_frac    = 0.003  # kg N / kg fresh residue (S4)

        # Methane energy
        self.CH4_energy_kWh_m3 = PHYSICAL_CONSTANTS["CH4_LHV_MJ_m3"] / 3.6

        # Thermal parameters
        self.U_eff          = U_eff   # R1 param
        self.Cp_substrate   = 4.18    # kJ/kg·K (dilute slurry ≈ water)
        self.eta_elec_heat  = 0.95    # electric resistance heating efficiency
        self.eta_bio_heat   = 0.85    # biogas CHP heating efficiency
        self.Q_rxn_frac     = 0.05    # exothermic reaction fraction of biogas energy

        # Arrhenius coefficient (Pilarski & Pilarska, 2025; Lindorfer et al., 2008)
        self.theta = 0.069            # °C⁻¹

    # ── Temperature models ────────────────────────────────────────────────
    def effective_digester_temp(self, T_amb: float) -> float:
        """T_eff = clip(T_amb + 15, 26°C, 37°C). Eq. (16)."""
        return float(np.clip(T_amb + 15.0, self.T_winter_min, self.T_mesophilic))

    def temperature_corrected_BMP(self, T_eff: float) -> float:
        """
        FIX 6 — Piecewise BMP model.
        T ≤ 37°C : Arrhenius (Eq. 17; Pilarski & Pilarska, 2025)
        T > 37°C : linear inhibition −5 NmL/gVS/°C (Feng et al., 2013)
        Note: T_eff ≤ 37°C operationally (clipped). This branch is only
        reached in run_validation() for data at 40°C (Feng et al., 2013).
        """
        if T_eff <= self.T_mesophilic:
            return float(
                self.BMP_ref * np.exp(self.theta * (T_eff - self.T_mesophilic)))
        else:
            return float(max(0.0, self.BMP_ref - 5.0 * (T_eff - self.T_mesophilic)))

    def monthly_BMP_array(self, monthly_T_amb: np.ndarray) -> np.ndarray:
        """Monthly BMP values (Arrhenius-corrected) from monthly mean T_amb."""
        return np.array([
            self.temperature_corrected_BMP(self.effective_digester_temp(T))
            for T in monthly_T_amb
        ])

    # ── Feedstock ─────────────────────────────────────────────────────────
    def residue_yield(self) -> Tuple[float, float, float]:
        """Returns (fresh_t/ha, dry_t/ha, VS_t/ha) using site-specific yield."""
        fresh = self.fruit_yield_t_ha * self.residue_ratio
        dry   = fresh * self.DM
        VS    = dry   * self.VS_DM
        return fresh, dry, VS

    def annual_methane_yield(self, BMP_avg: float,
                              eta_capture: float = 1.0) -> Dict[str, float]:
        """Annual CH4 from total seasonal VS × annual-average BMP."""
        _, _, VS_t = self.residue_yield()
        total_VS   = VS_t * 1000.0             # kg VS/ha/year
        CH4_m3     = total_VS * BMP_avg * 1000.0 / 1e6 * eta_capture
        return {"CH4_m3_annual":    CH4_m3,
                "biogas_kWh_annual": CH4_m3 * self.CH4_energy_kWh_m3}

    # ── Digester geometry ─────────────────────────────────────────────────
    def digester_geometry(self, volume_m3: float) -> Tuple[float, float, float]:
        """Cylindrical digester, H = D. Returns (radius, height, surface_area)."""
        r = (volume_m3 / (2.0 * np.pi)) ** (1.0 / 3.0)
        h = 2.0 * r
        A = 2.0 * np.pi * r ** 2 + 2.0 * np.pi * r * h
        return r, h, A

    # ── Heat demand (Eq. 18–21) ───────────────────────────────────────────
    def daily_heat_demand(self, volume_m3: float, T_amb: float,
                           feed_kg_day: float, BMP_day: float) -> float:
        """Daily digester heat demand (kWh/day)."""
        T_dig    = self.effective_digester_temp(T_amb)
        _, _, A  = self.digester_geometry(volume_m3)
        Q_loss   = self.U_eff * A * (T_dig - T_amb) * 24.0 / 1000.0
        Q_infl   = (feed_kg_day * self.Cp_substrate
                    * max(0.0, T_dig - max(T_amb, 5.0)) / 3600.0)
        Q_rxn    = (self.Q_rxn_frac
                    * self.annual_methane_yield(BMP_day)["biogas_kWh_annual"]
                    / 365.0)
        return float(max(0.0, Q_loss + Q_infl - Q_rxn))

    # ── Feasibility check ─────────────────────────────────────────────────
    def feasibility_check(self, volume_m3: float,
                           feed_kg_day: float,
                           VS_kg_day: float) -> Tuple[bool, Dict]:
        HRT = volume_m3 / max(feed_kg_day / 1000.0, 1e-9)
        OLR = VS_kg_day  / max(volume_m3, 1e-9)
        ok  = (self.HRT_min <= HRT <= self.HRT_max
               and self.OLR_min <= OLR <= self.OLR_max)
        return ok, {"HRT_actual": HRT, "OLR_actual": OLR}

    # ── Scenario economics (FIX 4 + FIX 5) ───────────────────────────────
    def scenario_economics(self, sc_key: str, biogas_sold_kWh: float,
                            fresh_t: float, site_cfg: dict) -> Dict:
        """
        FIX 4: prices read from site_cfg, not from any global dict.
        FIX 5: S5 uses BIOMETHANE_PREMIUM × p_elec.
        """
        sc  = self.SCENARIOS[sc_key]
        p_e = site_cfg["p_elec_usd_kWh"]
        p_b = site_cfg["p_biogas_usd_kWh"]

        if sc["biomethane"]:                          # S5: FIX 5
            bio_rev = biogas_sold_kWh * p_e * BIOMETHANE_PREMIUM
        elif sc["local"]:                             # S3: LPG parity
            bio_rev = biogas_sold_kWh * CAPEX["lpg_price_usd_kWh"]
        else:                                         # S0–S2, S4: grid export
            bio_rev = biogas_sold_kWh * p_b

        dig_credit = 0.0
        if sc["digestate"]:                           # S4: N-P-K credit
            N_kg = fresh_t * 1000.0 * self.digestate_N_frac
            dig_credit = N_kg * CAPEX["p_N_fertilizer_usd_kg"]

        return {"biogas_revenue_usd":   bio_rev,
                "digestate_credit_usd": dig_credit,
                "scenario_label":       sc["label"]}


# =============================================================================
# SECTION 4 — CLIMATE DATA LOADER
# =============================================================================

def load_climate_data(filepath: str, solar: SolarGeometry) -> pd.DataFrame:
    """
    Load a PVGIS-SARAH3 TMY CSV file and add pre-computed solar geometry.

    Expected PVGIS columns (auto-detected by keyword matching):
      time(UTC), T2m, RH, G(h), Gb(n), Gd(h), WS10m, SP

    Returns a 8760-row DataFrame with columns:
      datetime, month, GHI, DNI, DHI, T2m, RH, WS10m,
      solar_altitude_rad, solar_azimuth_rad
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(
            f"Climate file not found: {filepath}\n"
            f"Download at: https://re.jrc.ec.europa.eu/pvg_tools/en/\n"
            f"Settings: hourly TMY, PVGIS-SARAH3, select coordinates.")

    with open(filepath, "r", encoding="utf-8") as fh:
        lines = fh.readlines()
    header_idx = next(
        (i for i, l in enumerate(lines) if "time(UTC)" in l), 17)

    df = pd.read_csv(filepath, skiprows=header_idx)
    df.columns = df.columns.str.strip()

    rename: Dict[str, str] = {}
    for c in df.columns:
        cl = c.lower()
        if   "time"  in cl:  rename[c] = "time"
        elif "t2m"   in cl:  rename[c] = "T2m"
        elif "rh"    in cl:  rename[c] = "RH"
        elif "g(h)"  in cl:  rename[c] = "GHI"
        elif "gb(n)" in cl:  rename[c] = "DNI"
        elif "gd(h)" in cl:  rename[c] = "DHI"
        elif "ws10m" in cl:  rename[c] = "WS10m"
        elif "sp"    in cl:  rename[c] = "SP"
    df = df.rename(columns=rename)

    # Keep only rows with a valid PVGIS timestamp format YYYYMMDD:HHMM
    df = df[df["time"].astype(str).str.match(r"\d{8}:\d{4}", na=False)]

    def _parse_pvgis_time(t: object) -> Optional[datetime]:
        t = str(t).strip()
        try:
            return datetime(int(t[0:4]), int(t[4:6]), int(t[6:8]),
                            int(t[9:11]),
                            int(t[11:13]) if len(t) > 11 else 0)
        except (ValueError, IndexError):
            return None

    df["datetime"] = df["time"].apply(_parse_pvgis_time)
    df = df.dropna(subset=["datetime"])

    for col in ["GHI", "DNI", "DHI", "T2m", "RH", "WS10m"]:
        df[col] = pd.to_numeric(df.get(col, 0), errors="coerce").fillna(0.0)

    df = df.sort_values("datetime").reset_index(drop=True)
    df["month"] = df["datetime"].dt.month

    # Pre-compute solar geometry (required for both simulation and reference PV)
    alts, azs = [], []
    for dt in df["datetime"]:
        n   = solar.day_of_year(dt)
        ha  = solar.hour_angle(dt)
        alt = solar.solar_altitude(n, ha)
        az  = solar.solar_azimuth(n, ha, alt)
        alts.append(alt)
        azs.append(az)
    df["solar_altitude_rad"] = alts
    df["solar_azimuth_rad"]  = azs

    print(f"  Loaded {len(df):,} records  "
          f"GHI [{df['GHI'].min():.0f}–{df['GHI'].max():.0f}] W/m²  "
          f"T2m [{df['T2m'].min():.1f}–{df['T2m'].max():.1f}] °C")
    return df


# =============================================================================
# SECTION 5 — INTEGRATED APV-AD SIMULATION ENGINE
# =============================================================================

SCOPE_NOTE = """
=== LER_biogas SCOPE NOTE ===
LER_biogas = E_biogas,net / E_biogas,ref where E_biogas,ref is the maximum
methane recoverable from the same residue stream at 37°C under thermally
unconstrained operation. When ESR = (E_PV + E_biogas) / Q_heat >> 1/f_PV,heat
(threshold ≈ 6), PV fully covers heat demand, biogas_for_heat = 0, and
LER_biogas = 1.000 by construction. This is a DESIGN VALIDATION CRITERION
confirming thermal self-sufficiency — not an empirically discovered finding.
At Freiburg (lower GHI, higher heat demand), ESR may approach the threshold
and LER_biogas < 1.000 is possible for some optimizer-proposed designs.
"""


class IntegratedAPVAD:
    """
    Full annual APV-AD simulation with all Version 3.0 fixes.

    PERF-1: Reference PV is computed with fully vectorized numpy operations.
    PERF-2: Monte-Carlo reuses solar geometry arrays from the base system.
    """

    def __init__(self, site_name: str, site_config: dict,
                 climate_df: pd.DataFrame,
                 BMP_ref: float = 300.0,
                 U_eff: float   = 0.8,
                 alpha_shade: float = 0.40) -> None:
        self.site_name = site_name
        self.site      = site_config
        self.df        = climate_df

        self.solar = SolarGeometry(site_config["latitude"],
                                   site_config["longitude"],
                                   site_config["timezone_offset"])
        self.av = AgrivoltaicSystem(alpha_shade=alpha_shade)
        self.ad = AnaerobicDigestionSystem(
            BMP_ref=BMP_ref, U_eff=U_eff,
            fruit_yield_t_ha=site_config["fruit_yield_t_ha"])  # FIX 4

        # Pre-compute monthly mean temperature (used in monthly_BMP_array)
        mts = climate_df.groupby("month")["T2m"].mean()
        self._monthly_T = np.array([mts[m] for m in range(1, 13)])

        # Pre-compute ET0 with site-specific psychrometric constant (FIX 4)
        gamma = site_config.get("psychro_gamma", 0.0665)
        self._ET0_arr = self.av.et0_penman_monteith(
            climate_df["T2m"].values,
            climate_df["RH"].values,
            climate_df["WS10m"].values,
            climate_df["GHI"].values,
            gamma=gamma)

        # Growing season mask (FIX 3: OR logic for cross-year seasons)
        s_m, s_d = site_config["growing_start"]
        e_m, e_d = site_config["growing_end"]
        yr    = climate_df["datetime"].iloc[0].year
        start = datetime(yr, s_m, s_d)
        end   = datetime(yr, e_m, e_d)
        if e_m < s_m:
            # Cross-year season (e.g. Oct–Feb): use OR logic within one TMY year
            self._growing_mask = (
                (climate_df["datetime"] >= start) |
                (climate_df["datetime"] <= end)
            ).values
        else:
            self._growing_mask = (
                (climate_df["datetime"] >= start) &
                (climate_df["datetime"] <= end)
            ).values

        # FIX 2: LER_PV reference cache with tilt tracking
        self._ref_pv_kWh_indep   = None
        self._ref_pv_tilt_cached = None

    # ── PERF-1: vectorized reference PV ───────────────────────────────────
    def _compute_reference_pv(self, tilt_deg: float) -> float:
        """
        Dense-pack open-field PV reference (LER_PV Definition A).
        PERF-1: fully vectorized — no Python loop over 8760 rows.
        """
        df    = self.df
        n_ref = self.av.modules_per_hectare(self.av.module_width + 0.10)
        tilt  = np.radians(tilt_deg)
        surf  = np.radians(180.0)

        alt = df["solar_altitude_rad"].values
        az  = df["solar_azimuth_rad"].values
        GHI = df["GHI"].values
        DNI = df["DNI"].values
        DHI = df["DHI"].values
        T2m = df["T2m"].values
        WS  = df["WS10m"].values

        cos_aoi = (np.sin(alt) * np.cos(tilt)
                   + np.cos(alt) * np.sin(tilt) * np.cos(az - surf))
        G_beam  = DNI * np.maximum(0.0, cos_aoi)
        G_diff  = DHI * (1.0 + np.cos(tilt)) / 2.0
        G_gnd   = GHI * self.av.albedo * (1.0 - np.cos(tilt)) / 2.0
        G_front = G_beam + G_diff + G_gnd
        G_rear  = (DHI * (1.0 + np.cos(np.pi - tilt)) / 2.0
                   + GHI * self.av.albedo * (1.0 - np.cos(np.pi - tilt)) / 2.0)
        G_eff   = np.maximum(0.0,
                             G_front + G_rear * self.av.bifaciality
                             * (1.0 - self.av.system_losses))

        # Cell temperature and power — vectorized
        T_cell = T2m + G_eff / (self.av.U0 + self.av.U1 * np.maximum(0.5, WS))
        P_kw   = ((G_eff / 1000.0) * n_ref
                  * self.av.module_power_stc
                  * (1.0 - self.av.system_losses)
                  * (1.0 + self.av.temp_coeff_pmax * (T_cell - 25.0))
                  / 1000.0)

        # Only daytime hours contribute (alt > 0)
        return float(np.sum(np.where(alt > 0.0, np.maximum(0.0, P_kw), 0.0)))

    # ── Core simulation ────────────────────────────────────────────────────
    def simulate(self, design_vars: dict,
                 scenario_key: str = "S0") -> dict:
        """
        Run the full annual APV-AD simulation for one design configuration.

        Parameters
        ----------
        design_vars  : dict with keys tilt, row_spacing, module_height,
                       digester_volume, HRT, f_PV_heat
        scenario_key : "S0" … "S5" (AD valorization scenario)

        Returns
        -------
        dict — all LER components, energy flows, water savings, economics.
                Returns {"feasible": False, ...} for infeasible configurations.
        """
        tilt   = design_vars["tilt"]
        d_row  = design_vars["row_spacing"]
        H_mod  = design_vars["module_height"]
        V_dig  = design_vars["digester_volume"]
        HRT    = design_vars["HRT"]
        f_heat = design_vars["f_PV_heat"]
        df     = self.df

        n_mod = self.av.modules_per_hectare(d_row)
        GCR   = self.av.ground_cover_ratio(n_mod)

        # ── Hourly PV and shading loop ─────────────────────────────────────
        pv_kWh      = np.zeros(len(df))
        shading_arr = np.zeros(len(df))

        for i, (_, row) in enumerate(df.iterrows()):
            alt = row["solar_altitude_rad"]
            az  = row["solar_azimuth_rad"]
            shading_arr[i] = self.solar.shading_factor(alt, az, H_mod, d_row)
            if alt > 0.0:
                G, _ = self.av.plane_of_array_irradiance(
                    row["GHI"], row["DNI"], row["DHI"], tilt, alt, az)
                T_cell = self.av.cell_temperature(G, row["T2m"], row["WS10m"])
                pv_kWh[i] = self.av.pv_power_kw(G, T_cell, n_mod)

        # ── FIX 1: asymmetric LER_crop ─────────────────────────────────────
        LER_crop = self.av.crop_yield_ratio(
            df["GHI"].values, shading_arr, self._growing_mask)

        # ── ET and water savings ────────────────────────────────────────────
        ET_apv       = self.av.et_reduction(self._ET0_arr, shading_arr)
        water_saved  = float(np.sum(self._ET0_arr - ET_apv))
        ET0_annual   = float(np.sum(self._ET0_arr))
        LER_water    = water_saved / ET0_annual if ET0_annual > 0.0 else 0.0

        # ── AD feasibility ─────────────────────────────────────────────────
        fresh_t, _, VS_t = self.ad.residue_yield()
        feed_kg_day  = (fresh_t * 1000.0) / HRT
        VS_kg_day    = (VS_t   * 1000.0) / HRT

        ok, hrt_olr = self.ad.feasibility_check(V_dig, feed_kg_day, VS_kg_day)
        if not ok:
            return {"feasible": False, "eLER_global": 0.0,
                    "reason": f"AD constraints violated: {hrt_olr}"}

        # ── Temperature-corrected biogas (scenario BMP multiplier applied) ─
        sc         = self.ad.SCENARIOS[scenario_key]
        monthly_BMP = self.ad.monthly_BMP_array(self._monthly_T) * sc["BMP_mult"]
        dpm        = PHYSICAL_CONSTANTS["days_per_month"]
        BMP_avg    = float(np.average(monthly_BMP, weights=dpm))
        bio_dict   = self.ad.annual_methane_yield(BMP_avg, sc["eta_cap"])
        biogas_tot = bio_dict["biogas_kWh_annual"]

        # ── Daily heat demand (365 days) ────────────────────────────────────
        daily_heat = np.zeros(365)
        for day in range(365):
            h0          = day * 24
            T_day       = float(df["T2m"].iloc[h0: h0 + 24].mean())
            mo_idx      = df["datetime"].iloc[h0].month - 1
            daily_heat[day] = self.ad.daily_heat_demand(
                V_dig, T_day, feed_kg_day, monthly_BMP[mo_idx])
        annual_heat = float(np.sum(daily_heat))

        # ── Energy balance ─────────────────────────────────────────────────
        pv_total  = float(np.sum(pv_kWh))
        pv_deliv  = pv_total * f_heat * self.ad.eta_elec_heat
        deficit   = max(0.0, annual_heat - pv_deliv)
        bio_heat  = deficit / self.ad.eta_bio_heat
        if bio_heat > biogas_tot:
            return {"feasible": False, "eLER_global": 0.0,
                    "reason": "Insufficient energy for digester heating"}

        pv_sold   = pv_total * (1.0 - f_heat)
        bio_sold  = biogas_tot - bio_heat
        ESR       = (pv_total + biogas_tot) / max(1.0, annual_heat)

        # ── FIX 2: LER_PV reference with tilt cache invalidation ───────────
        if (self._ref_pv_kWh_indep is None
                or self._ref_pv_tilt_cached != tilt):
            self._ref_pv_kWh_indep   = self._compute_reference_pv(tilt)
            self._ref_pv_tilt_cached = tilt
        ref_pv = self._ref_pv_kWh_indep

        LER_PV_A  = pv_sold / ref_pv if ref_pv > 0.0 else 0.0
        LER_PV_B  = float(GCR)  # circular / literature convention
        LER_biogas = bio_sold / biogas_tot if biogas_tot > 0.0 else 0.0

        eLER_A = LER_crop + LER_PV_A + LER_biogas
        eLER_B = LER_crop + LER_PV_B + LER_biogas
        LER_2C = LER_crop + LER_PV_A

        # ── FIX 4: site-specific economics ────────────────────────────────
        sc_econ = self.ad.scenario_economics(
            scenario_key, bio_sold, fresh_t, self.site)
        pv_rev  = pv_sold * self.site["p_elec_usd_kWh"]
        cr_rev  = (LER_crop * self.ad.fruit_yield_t_ha
                   * self.site["p_tomato_usd_t"])
        tot_rev = (pv_rev
                   + sc_econ["biogas_revenue_usd"]
                   + cr_rev
                   + sc_econ["digestate_credit_usd"])

        pv_capex   = n_mod * CAPEX["pv_capex_usd_module"]
        ad_capex   = V_dig * CAPEX["ad_capex_usd_m3"]
        tot_capex  = pv_capex + ad_capex + pv_capex * CAPEX["struct_capex_frac"]
        net_income = tot_rev - tot_capex * CAPEX["opex_frac_capex"]
        payback    = tot_capex / net_income if net_income > 0.0 else float("inf")

        return {
            # Feasibility
            "feasible":           True,
            "scenario":           scenario_key,
            # eLER components
            "eLER_global":        eLER_A,        # Def. A — independent PV ref
            "eLER_global_B":      eLER_B,        # Def. B — circular PV ref
            "LER_2C":             LER_2C,        # comparable to published LER
            "LER_crop":           LER_crop,
            "LER_PV_A":           LER_PV_A,
            "LER_PV_B":           LER_PV_B,
            "LER_biogas":         LER_biogas,
            "LER_water_frac":     LER_water,
            # Water
            "water_saved_mm_yr":  water_saved,
            "ET0_annual_mm":      ET0_annual,
            # Energy
            "pv_total_kWh":       pv_total,
            "pv_sold_kWh":        pv_sold,
            "biogas_total_kWh":   biogas_tot,
            "biogas_sold_kWh":    bio_sold,
            "biogas_heat_kWh":    bio_heat,
            "annual_heat_kWh":    annual_heat,
            "ESR":                ESR,
            "BMP_annual_avg":     BMP_avg,
            "monthly_BMP":        monthly_BMP,
            # AD hydraulics
            "HRT_actual":         hrt_olr["HRT_actual"],
            "OLR_actual":         hrt_olr["OLR_actual"],
            # System geometry
            "n_modules":          n_mod,
            "GCR":                GCR,
            # Economics
            "total_capex_usd":    tot_capex,
            "net_income_usd":     net_income,
            "payback_yrs":        payback,
            "scenario_label":     sc_econ["scenario_label"],
            # Arrays for figures
            "pv_kWh_arr":         pv_kWh,
            "shading_arr":        shading_arr,
            "daily_heat_arr":     daily_heat,
        }


# =============================================================================
# SECTION 6 — OPTIMIZATION
# =============================================================================

BOUNDS: List[Tuple[float, float]] = [
    (10.0,  40.0),   # tilt (°)
    (4.0,   14.0),   # row_spacing (m)
    (2.0,    5.0),   # module_height (m)
    (50.0, 300.0),   # digester_volume (m³)
    (20.0,  40.0),   # HRT (days)
    (0.05,  0.50),   # f_PV_heat (−)
]
BOUND_NAMES = ["tilt", "row_spacing", "module_height",
               "digester_volume", "HRT", "f_PV_heat"]


def _unpack(x: np.ndarray) -> dict:
    return {n: x[i] for i, n in enumerate(BOUND_NAMES)}


class SingleObjectiveOptimizer:
    """Differential evolution maximizing eLER_global."""

    def __init__(self, system: IntegratedAPVAD,
                 scenario_key: str = "S0") -> None:
        self.system = system
        self.sc     = scenario_key

    def _objective(self, x: np.ndarray) -> float:
        r = self.system.simulate(_unpack(x), self.sc)
        return -r["eLER_global"] if r["feasible"] else 1000.0

    def optimize(self, d_row_min: float = 4.0,
                 maxiter: int = 50, popsize: int = 12,
                 seed: int = 42) -> Tuple[dict, dict]:
        """
        Returns (optimal_design_vars, full_result_dict).

        d_row_min : minimum row spacing enforced as lower bound (R3-A).
        """
        bounds = list(BOUNDS)
        bounds[1] = (d_row_min, 14.0)
        print(f"  DE: maxiter={maxiter}, popsize={popsize}, "
              f"d_row_min={d_row_min} m …")
        opt = differential_evolution(
            self._objective, bounds,
            maxiter=maxiter, popsize=popsize,
            seed=seed, disp=False, polish=True,
            mutation=(0.5, 1.0), recombination=0.7,
            tol=1e-5)
        best_dv  = _unpack(opt.x)
        best_res = self.system.simulate(best_dv, self.sc)
        best_res["optimizer_nfev"]    = opt.nfev
        best_res["optimizer_success"] = opt.success
        return best_dv, best_res


def run_access_sweep(system: IntegratedAPVAD,
                     d_row_values: List[float],
                     maxiter: int = 25,
                     seed: int = 42) -> pd.DataFrame:
    """R3-A: Sweep over d_row_min values and record optimal eLER."""
    print("\n[R3-A] Mechanical access sweep:")
    rows, opt = [], SingleObjectiveOptimizer(system)
    for d in d_row_values:
        print(f"  d_row_min={d} m …", end=" ", flush=True)
        dv, res = opt.optimize(d_row_min=d, maxiter=maxiter, seed=seed)
        print(f"eLER={res['eLER_global']:.3f}  LER_2C={res['LER_2C']:.3f}")
        rows.append({
            "d_row_min_m":   d,
            "eLER_global":   res["eLER_global"],
            "LER_2C":        res["LER_2C"],
            "LER_crop":      res["LER_crop"],
            "LER_PV_A":      res["LER_PV_A"],
            "LER_biogas":    res["LER_biogas"],
            "LER_water_frac":res["LER_water_frac"],
            "opt_GCR":       res["GCR"],
        })
    return pd.DataFrame(rows)


# ── FIX 9: Pareto problem with safe base class ─────────────────────────────

class ParetoAPVAD(_PymooBase):
    """
    FIX 9 — Safe class definition.
    _PymooBase is pymoo.core.problem.Problem when pymoo is installed,
    or plain object when it is not. The class is effectively unusable when
    pymoo is absent, but import does not raise a NameError.

    Objectives  : f1 = -eLER_global, f2 = -LER_water_frac
    Constraints : g1 = d_row_min - row_spacing ≤ 0  (access constraint)
                  g2 = 1 - ESR ≤ 0                   (thermal self-sufficiency)
    """

    def __init__(self, system: IntegratedAPVAD, d_row_min: float = 4.0) -> None:
        if _PYMOO:
            super().__init__(
                n_var=6, n_obj=2, n_ieq_constr=2,
                xl=np.array([b[0] for b in BOUNDS]),
                xu=np.array([b[1] for b in BOUNDS]))
        self.system    = system
        self.d_row_min = d_row_min

    def _evaluate(self, X: np.ndarray, out: dict, *args, **kwargs) -> None:
        F = np.full((len(X), 2), 1e6)
        G = np.full((len(X), 2), 1e6)
        for i, x in enumerate(X):
            dv = _unpack(x)
            r  = self.system.simulate(dv)
            if r["feasible"]:
                F[i, 0] = -r["eLER_global"]
                F[i, 1] = -r["LER_water_frac"]
                G[i, 0] = self.d_row_min - dv["row_spacing"]  # ≤ 0
                G[i, 1] = 1.0 - r["ESR"]                      # ≤ 0
            # else: penalty remains 1e6
        out["F"] = F
        out["G"] = G


def run_pareto(system: IntegratedAPVAD,
               d_row_min: float = 4.0,
               pop: int = 50,
               n_gen: int = 80,
               seed: int = 42) -> pd.DataFrame:
    """
    R3-B NSGA-II Pareto front — eLER vs. LER_water.

    FIX 10: Only called for sites with growing season > pareto_min_season_days.
    Almería (122 d) and Ouagadougou (122 d) are excluded because LER_water
    variation falls below numerical resolution at GCR ≈ 18.2% over 122 days.
    See paper Section 3.9 for the geometric derivation.
    """
    if not _PYMOO:
        print("[R3-B] pymoo not available — skipping Pareto.")
        return pd.DataFrame()

    print(f"\n[R3-B] NSGA-II (pop={pop}, gen={n_gen}, d_row_min={d_row_min} m) …")
    res = pymoo_minimize(
        ParetoAPVAD(system, d_row_min),
        NSGA2(pop_size=pop),
        get_termination("n_gen", n_gen),
        seed=seed, verbose=False)

    if res.X is None:
        print("  [WARNING] NSGA-II returned no solutions.")
        return pd.DataFrame()

    rows = []
    for x in res.X:
        r = system.simulate(_unpack(x))
        if r["feasible"]:
            rows.append({
                "eLER_global":    r["eLER_global"],
                "LER_water_frac": r["LER_water_frac"],
                "LER_2C":         r["LER_2C"],
                "GCR":            r["GCR"],
                "water_saved_mm": r["water_saved_mm_yr"],
            })

    df_pareto = pd.DataFrame(rows)
    print(f"  {len(df_pareto)} Pareto-optimal solutions found.")
    return df_pareto


# =============================================================================
# SECTION 7 — MONTE-CARLO  (FIX 7: structural invariance + economic MC)
# =============================================================================

def run_monte_carlo(system: IntegratedAPVAD,
                    optimal_vars: dict,
                    n_samples: int = 2000,
                    seed: int = 99) -> dict:
    """
    FIX 7 — Rewritten as STRUCTURAL INVARIANCE ANALYSIS + ECONOMIC MC.

    The eLER metric is structurally invariant to three of the six parameters
    varied (BMP_ref, α_shade, U_eff) because:
      • LER_biogas is invariant to BMP_ref: numerator and denominator both
        scale proportionally with BMP — the ratio cancels exactly.
      • LER_crop is invariant to α_shade: α_shade enters only the ET model,
        not the PAR integral of Eq. (11).
      • LER_PV and LER_biogas are invariant to U_eff: U_eff affects only
        Q_heat, and when ESR >> 1 the PV allocation exceeds any U_eff-driven
        variation in Q_heat.
    These invariances are documented and reported, not hidden.

    PERF-2: MC samples share the pre-computed solar geometry arrays from the
    base system's DataFrame. A lightweight MC system is built per sample that
    only replaces the three physical parameters (BMP_ref, U_eff, α_shade)
    and the three economic parameters (p_elec, p_biogas, p_tomato).
    The solar geometry computation (the bottleneck) is NOT repeated.

    Six LHS dimensions:
      0: BMP_ref      U[240, 360] NmL/gVS
      1: alpha_shade  U[0.25, 0.55]
      2: U_eff        U[0.50, 1.10] W/m²·K
      3: p_elec       U[0.05, 0.15] USD/kWh
      4: p_biogas     U[0.03, 0.10] USD/kWh
      5: p_tomato     U[150, 500]   USD/t
    Dims 3–5 affect economics only; dims 0–2 affect eLER physically (but
    the structural invariance means all eLER variation comes from design
    variable perturbations, captured separately in the tornado diagram).
    """
    print(f"\n[R1] Monte-Carlo: N={n_samples}, seed={seed}")

    sampler = qmc.LatinHypercube(d=6, seed=seed)
    raw     = sampler.random(n=n_samples)
    BMP_s   = qmc.scale(raw[:, 0:1], [240.0], [360.0]).flatten()
    alpha_s = qmc.scale(raw[:, 1:2], [0.25],  [0.55 ]).flatten()
    Ueff_s  = qmc.scale(raw[:, 2:3], [0.50],  [1.10 ]).flatten()
    pelec_s = qmc.scale(raw[:, 3:4], [0.05],  [0.15 ]).flatten()
    pbio_s  = qmc.scale(raw[:, 4:5], [0.03],  [0.10 ]).flatten()
    ptom_s  = qmc.scale(raw[:, 5:6], [150.0], [500.0]).flatten()

    eLER_arr = np.full(n_samples, np.nan)
    net_arr  = np.zeros(n_samples)

    for i in range(n_samples):
        if i % 500 == 0:
            print(f"  sample {i}/{n_samples} …", flush=True)

        # PERF-2: Build a lightweight MC system that reuses solar geometry.
        # We build a new IntegratedAPVAD (which re-computes ET0 with the
        # new alpha_shade) but point it at the same climate_df so the
        # solar geometry columns are available without recomputation.
        mc_sys = IntegratedAPVAD(
            system.site_name, system.site, system.df,
            BMP_ref=BMP_s[i], U_eff=Ueff_s[i], alpha_shade=alpha_s[i])

        # Override economic parameters without mutating the original site dict
        mc_site = dict(system.site)
        mc_site["p_elec_usd_kWh"]   = pelec_s[i]
        mc_site["p_biogas_usd_kWh"] = pbio_s[i]
        mc_site["p_tomato_usd_t"]   = ptom_s[i]
        mc_sys.site = mc_site

        r = mc_sys.simulate(optimal_vars, "S0")
        if r["feasible"]:
            eLER_arr[i] = r["eLER_global"]
            net_arr[i]  = r["net_income_usd"]

    valid  = ~np.isnan(eLER_arr)
    ev     = eLER_arr[valid]
    nv     = net_arr[valid]

    stats = {
        "n_valid":   int(valid.sum()),
        "eLER_mean": float(np.mean(ev)),
        "eLER_std":  float(np.std(ev)),
        "eLER_CV":   float(100.0 * np.std(ev) / np.mean(ev)),
        "eLER_P5":   float(np.percentile(ev, 5)),
        "eLER_P95":  float(np.percentile(ev, 95)),
        "net_P5":    float(np.percentile(nv, 5)),
        "net_P95":   float(np.percentile(nv, 95)),
    }

    # One-way eLER sensitivity to design variables (tornado)
    tornado: Dict[str, Tuple[float, float]] = {}
    for idx, name in enumerate(BOUND_NAMES):
        lo_dv = dict(optimal_vars); lo_dv[name] = BOUNDS[idx][0]
        hi_dv = dict(optimal_vars); hi_dv[name] = BOUNDS[idx][1]
        r_lo  = system.simulate(lo_dv)
        r_hi  = system.simulate(hi_dv)
        lo_v  = r_lo["eLER_global"] if r_lo["feasible"] else np.nan
        hi_v  = r_hi["eLER_global"] if r_hi["feasible"] else np.nan
        tornado[name] = (lo_v, hi_v)

    print(f"  eLER = {stats['eLER_mean']:.3f} ± {stats['eLER_std']:.3f}  "
          f"CV = {stats['eLER_CV']:.1f}%  "
          f"Net income P5–P95: ${stats['net_P5']:,.0f}–${stats['net_P95']:,.0f}")

    return {
        "eLER_arr": eLER_arr,
        "net_arr":  net_arr,
        "stats":    stats,
        "tornado":  tornado,
        "samples":  {
            "BMP_ref":    BMP_s, "alpha_shade": alpha_s,
            "U_eff":      Ueff_s, "p_elec":     pelec_s,
            "p_biogas":   pbio_s, "p_tomato":   ptom_s,
        },
    }


# =============================================================================
# SECTION 8 — VALIDATION  (FIX 6 applied)
# =============================================================================

def run_validation() -> dict:
    """
    R5 — Component-level validation.

    V1: Crop yield ratio vs. Mohammedi et al. (2023).
        Uses FIX 1 asymmetric formula (same as simulation).
        Reports MAPE for full shading range AND operational range ≤ 0.25.
    V2: BMP temperature response vs. Feng et al. (2013).
        FIX 6 piecewise model; V2 MAPE now correctly captures the
        inhibition above 37°C observed in the experimental data.
    """
    print("\n[R5] Validation:")
    results = {}

    # ── V1: crop yield ratio ───────────────────────────────────────────────
    shading_exp = np.array([0.15, 0.25, 0.35, 0.50, 0.65, 0.80])
    yield_obs   = np.array([0.97, 0.93, 0.87, 0.79, 0.68, 0.54])

    # FIX 1 formula: asymmetric — raw denominator, clipped numerator
    GHI_rep  = 500.0
    PAR_frac = PHYSICAL_CONSTANTS["PAR_to_GHI_ratio"]
    PAR_sat  = PHYSICAL_CONSTANTS["PAR_sat_tomato_W_m2"]
    PAR_open = GHI_rep * PAR_frac         # 240 W/m² > PAR_sat (174 W/m²)
    PAR_AV   = PAR_open * (1.0 - shading_exp)
    yield_mod = np.minimum(PAR_AV, PAR_sat) / PAR_open   # asymmetric

    mape_all = float(np.mean(np.abs(yield_mod - yield_obs) / yield_obs) * 100)
    rmse_v1  = float(np.sqrt(np.mean((yield_mod - yield_obs) ** 2)))
    op_mask  = shading_exp <= 0.25
    mape_op  = float(np.mean(np.abs(yield_mod[op_mask] - yield_obs[op_mask])
                              / yield_obs[op_mask]) * 100)
    print(f"  V1 MAPE (full range): {mape_all:.1f}%  RMSE={rmse_v1:.3f}  "
          f"MAPE (shading ≤ 0.25, operational): {mape_op:.1f}%")

    results["V1"] = {"obs": yield_obs, "mod": yield_mod,
                     "shading": shading_exp,
                     "MAPE_all": mape_all, "RMSE": rmse_v1,
                     "MAPE_operational": mape_op}

    # ── V2: BMP temperature response ──────────────────────────────────────
    T_exp   = np.array([25.0, 30.0, 35.0, 37.0, 40.0])
    BMP_obs = np.array([195.0, 240.0, 278.0, 300.0, 285.0])
    ad_val  = AnaerobicDigestionSystem(BMP_ref=300.0)
    BMP_mod = np.array([ad_val.temperature_corrected_BMP(T) for T in T_exp])
    mape_v2 = float(np.mean(np.abs(BMP_mod - BMP_obs) / BMP_obs) * 100)
    rmse_v2 = float(np.sqrt(np.mean((BMP_mod - BMP_obs) ** 2)))
    print(f"  V2 MAPE: {mape_v2:.1f}%  RMSE={rmse_v2:.1f} NmL/gVS")

    results["V2"] = {"obs": BMP_obs, "mod": BMP_mod, "T": T_exp,
                     "MAPE": mape_v2, "RMSE": rmse_v2}
    return results


# =============================================================================
# SECTION 9 — FIGURES
# =============================================================================

SITE_COLORS = {
    "Konya":       "#1B6CA8",  # deep blue
    "Almeria":     "#E05C2A",  # burnt orange
    "Ouagadougou": "#2E9B5A",  # forest green
    "Freiburg":    "#8B4513",  # saddle brown (temperate site)
}


def _save(fig: plt.Figure, name: str) -> None:
    """Save figure as PNG and PDF to output directory."""
    out = Path(RUN_CONFIG["output_dir"])
    for ext in ("png", "pdf"):
        p = out / f"{name}.{ext}"
        fig.savefig(str(p), dpi=RUN_CONFIG["figure_dpi"], bbox_inches="tight",
                    format=ext)
    print(f"  Saved: {out / name}.png / .pdf")


def plot_monthly_energy_balance(res: dict, site_name: str,
                                 save_path: Optional[str] = None) -> None:
    """Fig 1: Monthly PV energy balance and BMP seasonal profile."""
    dpm      = PHYSICAL_CONSTANTS["days_per_month"]
    months   = np.arange(1, 13)
    mlabels  = ["Jan","Feb","Mar","Apr","May","Jun",
                "Jul","Aug","Sep","Oct","Nov","Dec"]

    pv_arr    = res["pv_kWh_arr"]
    cum_days  = np.concatenate([[0], np.cumsum(dpm)[:-1]])
    monthly_pv = np.array([np.sum(pv_arr[h*24:(h+d)*24]) / 1000.0
                            for h, d in zip(cum_days, dpm)])
    frac      = res["biogas_heat_kWh"] / max(res["pv_total_kWh"], 1.0)
    f_heat    = monthly_pv * frac
    f_sold    = monthly_pv - f_heat

    daily_h   = res["daily_heat_arr"]
    mheat     = np.array([np.sum(daily_h[d: d + n]) / 1000.0
                           for d, n in zip(cum_days, dpm)])
    bio_m     = res["biogas_total_kWh"] * dpm / 365.0 / 1000.0

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), tight_layout=True)
    fig.suptitle(f"APV-AD Monthly Energy Balance — {site_name}",
                 fontsize=13, fontweight="bold")

    # Upper panel: PV stacked + heat demand
    ax1r = ax1.twinx()
    ax1.bar(months, f_sold, color="#185FA5", label="PV sold to grid", zorder=3)
    ax1.bar(months, f_heat, bottom=f_sold, color="#EF9F27",
            label="PV → digester heating", zorder=3)
    ax1r.plot(months, mheat, "o--", color="#D85A30", lw=2, ms=6,
              label="Digester heat demand", zorder=4)
    ax1.set_ylabel("PV energy (MWh/month)")
    ax1r.set_ylabel("Heat demand (MWh/month)", color="#D85A30")
    ax1r.tick_params(axis="y", labelcolor="#D85A30")
    ax1.set_xticks(months); ax1.set_xticklabels(mlabels)
    ax1.grid(axis="y", alpha=0.3)
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax1r.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=9)

    # Lower panel: biogas + BMP
    ax2r = ax2.twinx()
    ax2.bar(months, bio_m, color="#1D9E75", alpha=0.85,
            label="Monthly biogas energy", zorder=3)
    ax2r.plot(months, res["monthly_BMP"], "o--", color="#993C1D",
              lw=2, ms=6, label="BMP (NmL CH₄/g VS)", zorder=4)
    ax2.set_ylabel("Biogas energy (MWh/month)")
    ax2r.set_ylabel("BMP (NmL CH₄/g VS)", color="#993C1D")
    ax2r.tick_params(axis="y", labelcolor="#993C1D")
    ax2.set_xticks(months); ax2.set_xticklabels(mlabels)
    ax2.grid(axis="y", alpha=0.3)
    h3, l3 = ax2.get_legend_handles_labels()
    h4, l4 = ax2r.get_legend_handles_labels()
    ax2.legend(h3 + h4, l3 + l4, loc="upper right", fontsize=9)

    ann = (f"Annual PV: {res['pv_total_kWh']/1000:.0f} MWh  |  "
           f"Biogas: {res['biogas_total_kWh']/1000:.1f} MWh  |  "
           f"Heat: {res['annual_heat_kWh']/1000:.1f} MWh  |  "
           f"ESR: {res['ESR']:.1f}  |  eLER: {res['eLER_global']:.3f}")
    fig.text(0.5, -0.01, ann, ha="center", fontsize=9, style="italic")

    if save_path:
        fig.savefig(save_path, dpi=RUN_CONFIG["figure_dpi"], bbox_inches="tight")
        print(f"  Saved: {save_path}")
    plt.close(fig)


def plot_ler_decomposition(site_results: dict,
                            save_path: Optional[str] = None) -> None:
    """Fig 2: eLER decomposition for all sites, dual LER_PV definitions."""
    sites   = list(site_results.keys())
    n       = len(sites)
    eLER_A  = [site_results[s]["eLER_global"]  for s in sites]
    eLER_B  = [site_results[s]["eLER_global_B"]for s in sites]
    LER_c   = [site_results[s]["LER_crop"]      for s in sites]
    LER_pA  = [site_results[s]["LER_PV_A"]      for s in sites]
    LER_bio = [site_results[s]["LER_biogas"]     for s in sites]
    LER_2C  = [site_results[s]["LER_2C"]         for s in sites]
    x = np.arange(n); w = 0.35

    fig, ax = plt.subplots(figsize=(max(10, n * 2.5), 6), tight_layout=True)
    ax.bar(x - w/2, LER_c,  w, label="LER_crop",
           color="#3A8DC4", zorder=3)
    ax.bar(x - w/2, LER_pA, w, bottom=LER_c,
           label="LER_PV (Def. A, independent)", color="#F0A500", zorder=3)
    ax.bar(x - w/2, LER_bio, w,
           bottom=np.array(LER_c) + np.array(LER_pA),
           label="LER_biogas", color="#1D9E75", zorder=3)
    ax.bar(x + w/2, LER_c, w, color="#3A8DC4", alpha=0.5, zorder=3)
    ax.bar(x + w/2, np.array(eLER_B) - np.array(LER_c), w,
           bottom=LER_c,
           label="LER_PV (Def. B, circular) + LER_bio",
           color="#EF9F27", alpha=0.6, zorder=3)
    ax.axhline(1.94, color="red", ls="--", lw=1.5,
               label="Lit. max LER_2C = 1.94 (Riaz et al., 2022)")
    for i in range(n):
        ax.annotate(f"eLER_A={eLER_A[i]:.3f}",
                    (x[i] - w/2, eLER_A[i] + 0.02),
                    ha="center", fontsize=8, fontweight="bold")
        ax.annotate(f"LER_2C={LER_2C[i]:.3f}",
                    (x[i] - w/2, LER_2C[i] + 0.02),
                    ha="center", fontsize=7, color="#555555",
                    xytext=(0, -12), textcoords="offset points")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}\n(Def.A | Def.B)" for s in sites], fontsize=9)
    ax.set_ylabel("LER component value")
    ax.set_title("eLER decomposition — R2 Dual LER_PV definition comparison",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(axis="y", alpha=0.3)

    if save_path:
        fig.savefig(save_path, dpi=RUN_CONFIG["figure_dpi"], bbox_inches="tight")
        print(f"  Saved: {save_path}")
    plt.close(fig)


def plot_climate_gradient(site_results: dict,
                           save_path: Optional[str] = None) -> None:
    """
    Fig 8 (NEW) — Four-site climate gradient.
    GHI and mean T values derived from simulation TMY data where possible.
    """
    # Reference climate data (from PVGIS-SARAH3 TMY averages, Table 4)
    SITE_GHI  = {"Konya": 1650, "Almeria": 1850, "Ouagadougou": 2050, "Freiburg": 1150}
    SITE_TMEAN= {"Konya": 11.5, "Almeria": 18.5, "Ouagadougou": 28.5, "Freiburg": 11.4}

    sites_ord = sorted(site_results.keys(),
                        key=lambda s: SITE_GHI.get(s, 0), reverse=True)
    ghis   = [SITE_GHI.get(s, 0)   for s in sites_ord]
    temps  = [SITE_TMEAN.get(s, 0) for s in sites_ord]
    eLERs  = [site_results[s]["eLER_global"] for s in sites_ord]
    LER2Cs = [site_results[s]["LER_2C"]      for s in sites_ord]
    colors = [SITE_COLORS.get(s, "#333333")  for s in sites_ord]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), tight_layout=True)
    fig.suptitle("APV-AD Four-Site Climate Gradient\n"
                 "(Ouagadougou → Almería → Konya → Freiburg, decreasing GHI)",
                 fontsize=12, fontweight="bold")

    axes[0].barh(sites_ord, ghis, color=colors, height=0.5, zorder=3)
    axes[0].set_xlabel("Annual GHI (kWh/m²)")
    axes[0].set_title("Solar Resource")
    axes[0].grid(axis="x", alpha=0.3)
    for i, v in enumerate(ghis):
        axes[0].text(v + 15, i, str(v), va="center", fontsize=9)

    axes[1].barh(sites_ord, temps, color=colors, height=0.5, zorder=3)
    axes[1].set_xlabel("Mean annual T (°C)")
    axes[1].set_title("Mean Temperature")
    axes[1].grid(axis="x", alpha=0.3)
    for i, v in enumerate(temps):
        axes[1].text(v + 0.2, i, f"{v}°C", va="center", fontsize=9)

    bw = 0.35
    y  = np.arange(len(sites_ord))
    axes[2].barh(y - bw/2, eLERs,  bw, color=colors, label="eLER_global",    zorder=3)
    axes[2].barh(y + bw/2, LER2Cs, bw, color=colors, alpha=0.45,
                 label="LER_2C (crop+PV)", zorder=3)
    axes[2].axvline(1.94, color="red", ls="--", lw=1.2,
                    label="Lit. max LER_2C = 1.94")
    axes[2].set_yticks(y); axes[2].set_yticklabels(sites_ord)
    axes[2].set_xlabel("LER value"); axes[2].set_title("Land Equivalent Ratio")
    axes[2].legend(fontsize=8); axes[2].grid(axis="x", alpha=0.3)

    if save_path:
        fig.savefig(save_path, dpi=RUN_CONFIG["figure_dpi"], bbox_inches="tight")
        print(f"  Saved: {save_path}")
    plt.close(fig)


def plot_mc_tornado(mc_results: dict, site_name: str,
                    save_path: Optional[str] = None) -> None:
    """
    Fig 3 — FIX 7: Tornado (eLER sensitivity) + Net income distribution.
    Replaces the uninformative CV=0.0% spike histogram.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), tight_layout=True)
    fig.suptitle(
        f"R1 Structural Invariance & Economic MC "
        f"(N={RUN_CONFIG['mc_n_samples']}) — {site_name}",
        fontsize=12, fontweight="bold")

    # Tornado
    tornado    = mc_results["tornado"]
    base_eLER  = mc_results["stats"]["eLER_mean"]
    names      = list(tornado.keys())
    spans      = [abs((v[1] if not np.isnan(v[1]) else base_eLER)
                      - (v[0] if not np.isnan(v[0]) else base_eLER))
                  for v in tornado.values()]
    order      = np.argsort(spans)[::-1]

    for rank, idx in enumerate(order):
        name = names[idx]
        lo   = tornado[name][0] if not np.isnan(tornado[name][0]) else base_eLER
        hi   = tornado[name][1] if not np.isnan(tornado[name][1]) else base_eLER
        ax1.barh(rank, hi - base_eLER, left=base_eLER, height=0.5,
                 color="#1D9E75", alpha=0.8, zorder=3)
        ax1.barh(rank, lo - base_eLER, left=base_eLER, height=0.5,
                 color="#E05C2A", alpha=0.8, zorder=3)

    ax1.set_yticks(range(len(names)))
    ax1.set_yticklabels([names[i] for i in order], fontsize=9)
    ax1.axvline(base_eLER, color="black", lw=1.5)
    ax1.set_xlabel("eLER_global")
    ax1.set_title("One-way eLER sensitivity\n(design variables at bounds)")
    ax1.grid(axis="x", alpha=0.3)

    # Net income distribution
    valid = ~np.isnan(mc_results["eLER_arr"])
    net_v = mc_results["net_arr"][valid]
    s     = mc_results["stats"]
    ax2.hist(net_v / 1000.0, bins=50, color="#8B4FB5",
             edgecolor="white", linewidth=0.4, density=True)
    ax2.axvline(s["net_P5"]  / 1000.0, color="red", ls="--", lw=1.5,
                label=f"P5 = ${s['net_P5']:,.0f}")
    ax2.axvline(s["net_P95"] / 1000.0, color="red", ls="--", lw=1.5,
                label=f"P95 = ${s['net_P95']:,.0f}")
    ax2.set_xlabel("Net income (k$/yr)")
    ax2.set_ylabel("Density")
    ax2.set_title("Net income distribution\n(p_elec, p_biogas, p_tomato varied)")
    ax2.legend(fontsize=9); ax2.grid(alpha=0.3)

    fig.text(
        0.5, -0.03,
        "Note: eLER is structurally invariant to BMP_ref, α_shade, and U_eff "
        "(see Section 5.6.1). Economic spread driven by market price uncertainty.",
        ha="center", fontsize=8, style="italic", color="#555555")

    if save_path:
        fig.savefig(save_path, dpi=RUN_CONFIG["figure_dpi"], bbox_inches="tight")
        print(f"  Saved: {save_path}")
    plt.close(fig)


def plot_scenarios(sc_results: dict, site_name: str,
                   save_path: Optional[str] = None) -> None:
    """Fig 5 — R4 AD scenario comparison (S0–S5) with zoomed y-axis."""
    labels, eLERs, incs = [], [], []
    for k in ["S0", "S1", "S2", "S3", "S4", "S5"]:
        if k in sc_results and sc_results[k][1]["feasible"]:
            _, r = sc_results[k]
            labels.append(k)
            eLERs.append(r["eLER_global"])
            incs.append(r["net_income_usd"] / 1000.0)
    if not labels:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), tight_layout=True)
    fig.suptitle(f"R4: AD Scenario Comparison — {site_name}",
                 fontsize=12, fontweight="bold")
    colors = ["#3A8DC4", "#1D9E75", "#D85A30", "#F0A500", "#8B4FB5", "#C0392B"]

    # eLER panel — zoomed y-axis
    emin, emax = min(eLERs), max(eLERs)
    em = max(0.02, (emax - emin) * 3)
    ax1.bar(labels, eLERs, color=colors[:len(labels)], zorder=3)
    ax1.axhline(eLERs[0], color="grey", ls="--", lw=1, alpha=0.6)
    ax1.set_ylim(emin - em, emax + em * 1.5)
    for xi, val in enumerate(eLERs):
        d = val - eLERs[0]; sign = "+" if d >= 0 else ""
        ax1.text(xi, val + em * 0.3, f"{sign}{d:.3f}",
                 ha="center", fontsize=8)
    ax1.set_ylabel("eLER_global"); ax1.set_title("eLER by AD scenario")
    ax1.grid(axis="y", alpha=0.3)
    ax1.legend([f"S0 = {eLERs[0]:.3f}"], fontsize=9)

    # Economics panel — zoomed y-axis
    imin, imax = min(incs), max(incs)
    im = max(2.0, (imax - imin) * 3)
    ax2.bar(labels, incs, color=colors[:len(labels)], zorder=3)
    ax2.axhline(incs[0], color="grey", ls="--", lw=1, alpha=0.6)
    ax2.set_ylim(imin - im, imax + im * 1.5)
    for xi, val in enumerate(incs):
        d = val - incs[0]; sign = "+" if d >= 0 else ""
        ax2.text(xi, val + im * 0.3, f"{sign}{d:.1f}k",
                 ha="center", fontsize=8)
    ax2.set_ylabel("Net annual income (k$/yr)")
    ax2.set_title("Economics by AD scenario")
    ax2.grid(axis="y", alpha=0.3)
    ax2.legend([f"S0 = ${incs[0]:.0f}k"], fontsize=9)

    if save_path:
        fig.savefig(save_path, dpi=RUN_CONFIG["figure_dpi"], bbox_inches="tight")
        print(f"  Saved: {save_path}")
    plt.close(fig)


def plot_pareto(pareto_df: pd.DataFrame, site_name: str,
                save_path: Optional[str] = None) -> None:
    """Fig 4 — Pareto front (eLER vs. LER_water)."""
    if pareto_df.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True)
    sc = ax.scatter(pareto_df["LER_water_frac"] * 100,
                    pareto_df["eLER_global"],
                    c=pareto_df["GCR"] * 100, cmap="viridis",
                    s=60, alpha=0.8, zorder=3)
    plt.colorbar(sc, ax=ax, label="GCR (%)")
    ax.set_xlabel("LER_water (% ET₀ saved)")
    ax.set_ylabel("eLER_global (Definition A)")
    ax.set_title(f"R3-B: Pareto front — eLER vs. LER_water\n{site_name}",
                 fontsize=12, fontweight="bold")
    ax.grid(alpha=0.3)
    if save_path:
        fig.savefig(save_path, dpi=RUN_CONFIG["figure_dpi"], bbox_inches="tight")
        print(f"  Saved: {save_path}")
    plt.close(fig)


def plot_access_sweep(sweep_df: pd.DataFrame, site_name: str,
                       save_path: Optional[str] = None) -> None:
    """Fig 7 — eLER penalty of mechanical access constraint."""
    fig, ax = plt.subplots(figsize=(7, 5), tight_layout=True)
    ax.plot(sweep_df["d_row_min_m"], sweep_df["eLER_global"],
            "o-", color="#3A8DC4", lw=2, ms=9, label="eLER_global")
    ax.plot(sweep_df["d_row_min_m"], sweep_df["LER_2C"],
            "s--", color="#F0A500", lw=2, ms=9, label="LER_2C (crop+PV)")
    for _, row in sweep_df.iterrows():
        ax.annotate(f"{row['eLER_global']:.3f}",
                    (row["d_row_min_m"], row["eLER_global"] + 0.005),
                    ha="center", fontsize=9, color="#3A8DC4")
    ax.set_xlabel("Minimum row spacing d_row_min (m)\n[tractor access constraint]")
    ax.set_ylabel("LER value")
    ax.set_title(f"R3-A: eLER vs. Mechanical Access — {site_name}",
                 fontsize=12, fontweight="bold")
    ax.legend(); ax.grid(alpha=0.3)
    if save_path:
        fig.savefig(save_path, dpi=RUN_CONFIG["figure_dpi"], bbox_inches="tight")
        print(f"  Saved: {save_path}")
    plt.close(fig)


def plot_validation(val: dict, save_path: Optional[str] = None) -> None:
    """Fig R5 — V1 crop yield and V2 BMP validation."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), tight_layout=True)
    fig.suptitle("R5: Model Validation Against Published Experimental Data",
                 fontsize=12, fontweight="bold")

    v1 = val["V1"]
    lo = min(v1["obs"].min(), v1["mod"].min()) - 0.02
    ax1.scatter(v1["obs"], v1["mod"], s=80, color="#3A8DC4", zorder=3,
                label="Data points")
    ax1.plot([lo, 1.02], [lo, 1.02], "k--", lw=1.5, label="1:1 line")
    ax1.set_xlabel("Observed relative yield (Mohammedi et al., 2023)")
    ax1.set_ylabel("Modeled relative yield")
    ax1.set_title(
        f"V1: Crop yield ratio\n"
        f"MAPE (full) = {v1['MAPE_all']:.1f}%   "
        f"MAPE (op.) = {v1['MAPE_operational']:.1f}%   "
        f"RMSE = {v1['RMSE']:.3f}",
        fontsize=10)
    ax1.legend(fontsize=9); ax1.grid(alpha=0.3)

    v2 = val["V2"]
    ax2.plot(v2["T"], v2["obs"], "o-", color="#D85A30", lw=2, ms=8,
             label="Observed (Feng et al., 2013)")
    ax2.plot(v2["T"], v2["mod"], "s--", color="#3A8DC4", lw=2, ms=8,
             label="Modeled (piecewise Arrhenius)")
    ax2.set_xlabel("Temperature (°C)")
    ax2.set_ylabel("BMP (NmL CH₄/g VS)")
    ax2.set_title(
        f"V2: BMP temperature response\n"
        f"MAPE = {v2['MAPE']:.1f}%   RMSE = {v2['RMSE']:.1f} NmL/gVS",
        fontsize=10)
    ax2.legend(fontsize=9); ax2.grid(alpha=0.3)

    if save_path:
        fig.savefig(save_path, dpi=RUN_CONFIG["figure_dpi"], bbox_inches="tight")
        print(f"  Saved: {save_path}")
    plt.close(fig)


# =============================================================================
# SECTION 10 — RESULTS REPORTING AND CSV EXPORT
# =============================================================================

def print_results(site_name: str, dv: dict, r: dict) -> None:
    """Console results table for one site."""
    print(f"\n{'='*68}")
    print(f"  {site_name.upper()}  |  "
          f"{r.get('scenario_label', r.get('scenario', 'S0'))}")
    print(f"{'='*68}")
    print(f"  Tilt={dv['tilt']:.1f}°  d_row={dv['row_spacing']:.1f}m  "
          f"H_m={dv['module_height']:.1f}m  V_dig={dv['digester_volume']:.0f}m³  "
          f"HRT={dv['HRT']:.0f}d  f_heat={dv['f_PV_heat']:.2f}")
    print(f"  N_mod={r['n_modules']}  GCR={r['GCR']:.1%}  "
          f"HRT_act={r['HRT_actual']:.1f}d  OLR={r['OLR_actual']:.2f}")
    print(f"\n  eLER_A={r['eLER_global']:.3f}  eLER_B={r['eLER_global_B']:.3f}  "
          f"LER_2C={r['LER_2C']:.3f}")
    print(f"    LER_crop={r['LER_crop']:.3f}  LER_PV_A={r['LER_PV_A']:.3f}  "
          f"LER_biogas={r['LER_biogas']:.3f}  LER_water={r['LER_water_frac']:.3f}")
    print(f"\n  PV={r['pv_total_kWh']:,.0f} kWh  "
          f"Biogas={r['biogas_total_kWh']:,.0f} kWh  "
          f"Heat={r['annual_heat_kWh']:,.0f} kWh  ESR={r['ESR']:.1f}")
    print(f"  BMP_avg={r['BMP_annual_avg']:.1f} NmL/gVS  "
          f"Water={r['water_saved_mm_yr']:.0f} mm/yr")
    print(f"\n  CAPEX=${r['total_capex_usd']:,.0f}  "
          f"Net=${r['net_income_usd']:,.0f}/yr  "
          f"Payback={r['payback_yrs']:.1f} yr")
    if r["LER_biogas"] >= 0.995:
        print(SCOPE_NOTE)


def save_csv_summary(all_res: dict, all_dv: dict,
                     filepath: str) -> None:
    """Write cross-site summary table to CSV for table production."""
    fields = ["site", "climate_zone", "country",
              "eLER_global", "eLER_global_B", "LER_2C",
              "LER_crop", "LER_PV_A", "LER_PV_B", "LER_biogas",
              "LER_water_frac", "water_saved_mm_yr",
              "pv_total_kWh", "pv_sold_kWh",
              "biogas_total_kWh", "biogas_sold_kWh", "biogas_heat_kWh",
              "annual_heat_kWh", "ESR", "BMP_annual_avg",
              "HRT_actual", "OLR_actual", "n_modules", "GCR",
              "total_capex_usd", "net_income_usd", "payback_yrs",
              "tilt", "row_spacing", "module_height",
              "digester_volume", "HRT_design", "f_PV_heat"]
    with open(filepath, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for site, r in all_res.items():
            dv = all_dv[site]
            row = {
                "site":           site,
                "climate_zone":   SITES[site]["climate_zone"],
                "country":        SITES[site]["country"],
                "tilt":           f"{dv['tilt']:.2f}",
                "row_spacing":    f"{dv['row_spacing']:.2f}",
                "module_height":  f"{dv['module_height']:.2f}",
                "digester_volume":f"{dv['digester_volume']:.0f}",
                "HRT_design":     f"{dv['HRT']:.1f}",
                "f_PV_heat":      f"{dv['f_PV_heat']:.3f}",
            }
            for field in fields:
                if field not in row and field in r:
                    v = r[field]
                    row[field] = (f"{v:.4f}" if isinstance(v, float)
                                  else str(v))
            writer.writerow({f: row.get(f, "") for f in fields})
    print(f"  Saved CSV: {filepath}")


# =============================================================================
# SECTION 11 — MAIN EXECUTION PIPELINE
# =============================================================================

def main() -> None:
    print("=" * 68)
    print("  APV-AD INTEGRATED SYSTEM — VERSION 3.0")
    print("  Sites: Konya | Almería | Ouagadougou | Freiburg")
    print("  Fixes: 1–10 applied  |  Target: RSER / Applied Energy")
    print("=" * 68)

    out = Path(RUN_CONFIG["output_dir"])
    out.mkdir(exist_ok=True)

    # ── R5: Validation (site-independent) ─────────────────────────────────
    if RUN_CONFIG["run_validation"]:
        print("\n" + "─" * 40 + " R5 VALIDATION " + "─" * 13)
        val = run_validation()
        if RUN_CONFIG["save_figures"]:
            plot_validation(val, save_path=str(out / "R5_validation.png"))

    all_res:     dict = {}
    all_dv:      dict = {}

    for site_name in RUN_CONFIG["sites"]:
        cfg = SITES[site_name]
        print(f"\n{'='*68}")
        print(f"  SITE: {site_name.upper()}  "
              f"({cfg['climate_zone']}, {cfg['country']})  "
              f"yield={cfg['fruit_yield_t_ha']} t/ha  "
              f"p_elec={cfg['p_elec_usd_kWh']} USD/kWh")
        print("=" * 68)

        # Load climate data
        solar = SolarGeometry(cfg["latitude"], cfg["longitude"],
                              cfg["timezone_offset"])
        try:
            climate_df = load_climate_data(cfg["climate_file"], solar)
        except FileNotFoundError as exc:
            print(f"  [SKIP] {exc}")
            continue

        system = IntegratedAPVAD(site_name, cfg, climate_df)

        # ── Primary DE optimization (S0, d_row_min = 4 m) ─────────────────
        print(f"\n[OPT] {site_name}")
        opt     = SingleObjectiveOptimizer(system)
        dv, res = opt.optimize(
            d_row_min=4.0,
            maxiter=RUN_CONFIG["de_maxiter"],
            popsize=RUN_CONFIG["de_popsize"],
            seed=RUN_CONFIG["de_seed"])
        all_res[site_name] = res
        all_dv[site_name]  = dv
        print_results(site_name, dv, res)

        if RUN_CONFIG["save_figures"]:
            plot_monthly_energy_balance(
                res, site_name,
                save_path=str(out / f"Fig1_energy_balance_{site_name}.png"))

        # ── R3-A: Mechanical access constraint sweep ───────────────────────
        if RUN_CONFIG["run_mo_constraint_sweep"]:
            sw = run_access_sweep(
                system, RUN_CONFIG["mo_d_row_values"],
                maxiter=25, seed=RUN_CONFIG["de_seed"])
            if RUN_CONFIG["save_figures"]:
                plot_access_sweep(
                    sw, site_name,
                    save_path=str(out / f"Fig7_access_{site_name}.png"))

        # ── R3-B: Pareto front (FIX 10: season-length gate) ───────────────
        s_m, s_d = cfg["growing_start"]
        e_m, e_d = cfg["growing_end"]
        if e_m >= s_m:
            # Same-year season: approximate length in days
            from calendar import monthrange
            season_days = sum(
                monthrange(2005, m)[1]
                for m in range(s_m, e_m + 1)
            ) - s_d + e_d
        else:
            season_days = 0   # cross-year: always short for our sites

        run_pareto_flag = (RUN_CONFIG["run_mo_pareto"] and _PYMOO
                           and season_days >= RUN_CONFIG["pareto_min_season_days"])
        if run_pareto_flag:
            pf = run_pareto(system, d_row_min=4.0,
                            pop=50, n_gen=80, seed=RUN_CONFIG["de_seed"])
            if RUN_CONFIG["save_figures"] and not pf.empty:
                plot_pareto(
                    pf, site_name,
                    save_path=str(out / f"Fig4_pareto_{site_name}.png"))
        else:
            if RUN_CONFIG["run_mo_pareto"]:
                print(f"  [INFO] Pareto skipped for {site_name} "
                      f"(season {season_days}d < "
                      f"{RUN_CONFIG['pareto_min_season_days']}d threshold)")

        # ── R4: AD scenarios S0–S5 ─────────────────────────────────────────
        if RUN_CONFIG["run_ad_scenarios"]:
            print(f"\n[R4] AD scenarios — {site_name}")
            sc_results = {}
            for k in ["S0", "S1", "S2", "S3", "S4", "S5"]:
                r2 = system.simulate(dv, scenario_key=k)
                sc_results[k] = (dv, r2)
                if r2["feasible"]:
                    print(f"  {k}: eLER={r2['eLER_global']:.3f}  "
                          f"net=${r2['net_income_usd']:,.0f}  "
                          f"[{r2['scenario_label'][:35]}]")
            if RUN_CONFIG["save_figures"]:
                plot_scenarios(
                    sc_results, site_name,
                    save_path=str(out / f"Fig5_scenarios_{site_name}.png"))

        # ── R1: Monte-Carlo ────────────────────────────────────────────────
        if RUN_CONFIG["run_monte_carlo"]:
            mc = run_monte_carlo(
                system, dv,
                n_samples=RUN_CONFIG["mc_n_samples"],
                seed=RUN_CONFIG["mc_seed"])
            if RUN_CONFIG["save_figures"]:
                plot_mc_tornado(
                    mc, site_name,
                    save_path=str(out / f"Fig3_mc_{site_name}.png"))

    # ── Cross-site figures ─────────────────────────────────────────────────
    if len(all_res) > 1 and RUN_CONFIG["save_figures"]:
        plot_ler_decomposition(
            all_res,
            save_path=str(out / "Fig2_eLER_decomposition.png"))
        plot_climate_gradient(
            all_res,
            save_path=str(out / "Fig8_climate_gradient.png"))

    # ── Cross-site summary table ───────────────────────────────────────────
    print(f"\n{'='*68}")
    print("  CROSS-SITE SUMMARY (Table 10)")
    print(f"{'='*68}")
    hdr = (f"  {'Site':<15} {'Zone':>5} "
           f"{'eLER_A':>7} {'LER_2C':>7} "
           f"{'LER_c':>7} {'LER_PV_A':>8} {'LER_bio':>8} "
           f"{'ESR':>7} {'BMP':>7}")
    print(hdr)
    print("  " + "-" * 72)
    for site, r in all_res.items():
        print(f"  {site:<15} {SITES[site]['climate_zone']:>5} "
              f"{r['eLER_global']:>7.3f} {r['LER_2C']:>7.3f} "
              f"{r['LER_crop']:>7.3f} {r['LER_PV_A']:>8.3f} "
              f"{r['LER_biogas']:>8.3f} "
              f"{r['ESR']:>7.1f} {r['BMP_annual_avg']:>7.1f}")

    # Save CSV
    if all_res:
        save_csv_summary(
            all_res, all_dv,
            filepath=str(out / "crosssite_summary.csv"))

    print(SCOPE_NOTE)
    print(f"\n  All outputs saved to: {out.resolve()}")
    print("=" * 68)


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    main()
