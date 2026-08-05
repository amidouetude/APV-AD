# =============================================================================
# 00_config.py — Central Configuration for APV-AD Optimization Project
# =============================================================================
# All physical constants, site parameters, model settings, file paths,
# and 4E analysis factors are defined here.
# Every other notebook imports this file at the top:
#     from config_00 import SITES, PARAMS, PATHS, DE_SETTINGS, FACTORS_4E
#
# DO NOT hardcode any constant in other notebooks.
# If a value changes, change it here ONLY — it propagates everywhere.
# =============================================================================

from pathlib import Path

# =============================================================================
# 0. PROJECT PATHS
# =============================================================================

ROOT        = Path(__file__).resolve().parent.parent   # APV_AD_Project/
DATA_DIR    = ROOT / "data"
OUTPUTS_DIR = ROOT / "outputs"
CSV_DIR     = OUTPUTS_DIR / "csv"
FIG_DIR     = OUTPUTS_DIR / "figures"
PAPER_DIR   = ROOT / "paper" / "figures"

# Create output directories if they do not exist
for _dir in [CSV_DIR, FIG_DIR, PAPER_DIR]:
    _dir.mkdir(parents=True, exist_ok=True)

# =============================================================================
# 1. PVGIS TMY DATA FILE NAMES
#    Exact filenames as downloaded from PVGIS-SARAH3 and placed in data/
# =============================================================================

PVGIS_FILES = {
    "Konya":        DATA_DIR / "konya_climate.csv",
    "Almeria":      DATA_DIR / "alemeria_climate.csv",   # original spelling preserved
    "Ouagadougou":  DATA_DIR / "ouagadougou_climate.csv",
    "Freiburg":     DATA_DIR / "freiburg_climate.csv",
}

# =============================================================================
# 2. SITE PARAMETERS
#    One dictionary per site. All physical, agronomic, and economic
#    parameters that vary by location are stored here.
# =============================================================================

SITES = {

    # ── Konya, Turkey (BSk — cold semi-arid) ─────────────────────────────────
    "Konya": {
        # Geography
        "lat":              37.87,          # degrees North
        "lon":              32.49,          # degrees East
        "elevation_m":      1016,           # m above sea level
        "koppen":           "BSk",
        "timezone_utc":     3,              # UTC+3

        # Growing season — tomato (Solanum lycopersicum)
        "planting_date":    "04-25",        # MM-DD
        "harvest_date":     "10-30",        # MM-DD
        "season_days":      189,

        # Agronomic parameters
        "fruit_yield_t_ha": 60.0,           # fresh fruit yield (t/ha)
        "R_res":            0.80,           # residue-to-fruit ratio (Mohammedi 2023)
        "DM":               0.06,           # dry matter fraction — tomato residues (Mohammedi 2023)
        "VS_fraction":      0.85,           # volatile solids fraction of DM
        "PAR_sat":          174.0,          # W/m² — tomato light saturation point
        # Manure co-substrate (cattle, ~2 heads/ha baseline)
        "manure_input_t_ha":  10.0,      # t fresh manure/ha/yr (2-3 cattle)
        "manure_DM":           0.20,      # dry matter fraction — cattle manure
        "manure_VS_fraction":  0.80,      # VS/DM — cattle manure (Pilarski 2025)

        # Economic parameters
        "p_tomato_usd_t":   300.0,          # tomato farm gate price (USD/t)
        "p_elec_usd_kWh":   0.06,           # electricity sell price (USD/kWh)
        "p_biogas_usd_kWh": 0.08,           # baseline biogas price (USD/kWh)
        "discount_rate":    0.08,           # annual discount rate (8%)

        # 4E — Environmental factors
        "grid_ef_tCO2_MWh": 0.438,         # grid emission factor (IEA 2023 — Turkey)
        "carbon_price_usd_tCO2": 15.0,     # voluntary carbon market price (USD/tCO2)
    },

    # ── Almería, Spain (BSh — hot semi-arid) ─────────────────────────────────
    "Almeria": {
        # Geography
        "lat":              36.83,
        "lon":              -2.46,
        "elevation_m":      20,
        "koppen":           "BSh",
        "timezone_utc":     1,              # UTC+1 (CET)

        # Growing season
        "planting_date":    "03-01",
        "harvest_date":     "06-30",
        "season_days":      122,

        # Agronomic parameters
        "fruit_yield_t_ha": 60.0,
        "R_res":            0.80,
        "DM":               0.06,           # dry matter fraction — tomato residues (Mohammedi 2023)
        "VS_fraction":      0.85,
        "PAR_sat":          174.0,
        # Manure co-substrate (cattle, ~2 heads/ha baseline)
        "manure_input_t_ha":  10.0,      # t fresh manure/ha/yr (2-3 cattle)
        "manure_DM":           0.20,      # dry matter fraction — cattle manure
        "manure_VS_fraction":  0.80,      # VS/DM — cattle manure (Pilarski 2025)

        # Economic parameters
        "p_tomato_usd_t":   400.0,
        "p_elec_usd_kWh":   0.09,
        "p_biogas_usd_kWh": 0.09,
        "discount_rate":    0.06,           # 6%

        # 4E — Environmental factors
        "grid_ef_tCO2_MWh": 0.167,         # Red Electrica de Espana 2023
        "carbon_price_usd_tCO2": 65.0,     # EU ETS 2024 average (EUR ≈ USD here)
    },

    # ── Ouagadougou, Burkina Faso (BSh — tropical semi-arid) ─────────────────
    "Ouagadougou": {
        # Geography
        "lat":              12.37,
        "lon":              -1.53,
        "elevation_m":      303,
        "koppen":           "BSh",
        "timezone_utc":     0,              # UTC+0

        # Growing season
        "planting_date":    "07-01",
        "harvest_date":     "10-30",
        "season_days":      122,

        # Agronomic parameters
        "fruit_yield_t_ha": 60.0,
        "R_res":            0.80,
        "DM":               0.06,           # dry matter fraction — tomato residues (Mohammedi 2023)
        "VS_fraction":      0.85,
        "PAR_sat":          174.0,
        # Manure co-substrate (cattle, ~2 heads/ha baseline)
        "manure_input_t_ha":  10.0,      # t fresh manure/ha/yr (2-3 cattle)
        "manure_DM":           0.20,      # dry matter fraction — cattle manure
        "manure_VS_fraction":  0.80,      # VS/DM — cattle manure (Pilarski 2025)

        # Economic parameters
        "p_tomato_usd_t":   250.0,
        "p_elec_usd_kWh":   0.12,
        "p_biogas_usd_kWh": 0.08,
        "discount_rate":    0.10,           # 10%

        # 4E — Environmental factors
        "grid_ef_tCO2_MWh": 0.632,         # IPCC default Sub-Saharan diesel grid
        "carbon_price_usd_tCO2": 12.0,     # Gold Standard voluntary market
    },

    # ── Freiburg im Breisgau, Germany (Cfb — oceanic temperate) ──────────────
    "Freiburg": {
        # Geography
        "lat":              47.99,
        "lon":              7.85,
        "elevation_m":      278,
        "koppen":           "Cfb",
        "timezone_utc":     1,              # UTC+1 (CET)

        # Growing season
        "planting_date":    "05-01",
        "harvest_date":     "10-15",
        "season_days":      167,

        # Agronomic parameters
        "fruit_yield_t_ha": 30.0,           # German open-field statistics (DLG 2023)
        "R_res":            0.80,
        "DM":               0.06,           # dry matter fraction — tomato residues (Mohammedi 2023)
        "VS_fraction":      0.85,
        "PAR_sat":          174.0,
        # Manure co-substrate (cattle, ~2 heads/ha baseline)
        "manure_input_t_ha":  10.0,      # t fresh manure/ha/yr (2-3 cattle)
        "manure_DM":           0.20,      # dry matter fraction — cattle manure
        "manure_VS_fraction":  0.80,      # VS/DM — cattle manure (Pilarski 2025)

        # Economic parameters
        "p_tomato_usd_t":   350.0,
        "p_elec_usd_kWh":   0.09,
        "p_biogas_usd_kWh": 0.09,
        "discount_rate":    0.05,           # 5%

        # 4E — Environmental factors
        "grid_ef_tCO2_MWh": 0.380,         # Umweltbundesamt 2023
        "carbon_price_usd_tCO2": 65.0,     # EU ETS 2024
    },
}

# =============================================================================
# 3. PHYSICAL MODEL PARAMETERS
#    Fixed constants shared across all sites.
#    Source references match the equations in amidoumaiga_seminar2.tex
# =============================================================================

PARAMS = {

    # ── PV module — Jinko Tiger Neo 550Wp bifacial ────────────────────────────
    "P_STC_Wp":         550.0,      # rated power at STC (W)
    "bifaciality":      0.70,       # rear-to-front power ratio (phi)
    "beta_p":          -0.0035,     # temperature coefficient (K⁻¹)
    "eta_loss":         0.12,       # total system losses (wiring, soiling, mismatch)
    "rho_ground":       0.25,       # ground albedo for rear irradiance

    # ── Faiman thermal model (Faiman 2008) ────────────────────────────────────
    "U0":               25.0,       # W/m²·K — free convection coefficient
    "U1":               6.84,       # W·s/m³·K — wind-dependent coefficient

    # ── Crop yield model (McCree 1971; Mohammedi 2023) ────────────────────────
    "f_PAR":            0.48,       # fraction of GHI that is PAR
    "alpha_shade":      0.40,       # ET reduction coefficient under shade

    # ── Anaerobic digestion model ─────────────────────────────────────────────
    "BMP_ref":          300.0,      # NmL CH4/g VS at T_ref = 37°C (Lallement 2023)
    "T_ref_degC":       37.0,       # reference digester temperature (°C)
    "theta_arrhenius":  0.069,      # Arrhenius coefficient (°C⁻¹) (Pilarski 2025)
    "T_dig_min_degC":   26.0,       # minimum effective digester temperature (°C)
    "T_amb_boost_degC": 5.0,        # thermal boost above ambient (°C)
    "U_eff_W_m2K":      0.8,        # digester wall heat transfer coefficient
    "CH4_fraction":     0.60,       # methane fraction in biogas (vol/vol)
    "LHV_CH4_MJ_m3":   35.8,       # lower heating value of methane (MJ/Nm³)
    "eta_elec_heat":    0.95,       # electrical resistance heating efficiency
    "rho_substrate":    1020.0,     # substrate density (kg/m³)
    "Cp_substrate":     4.2,        # substrate specific heat (kJ/kg·K)
    "T_inlet_degC":     15.0,       # substrate inlet temperature (°C)
    "Q_reaction_frac":  0.05,       # exothermic reaction heat as fraction of E_biogas

    # ── 4E analysis — universal constants ────────────────────────────────────
    "ef_natural_gas_tCO2_MWh": 0.202,   # natural gas emission factor (IPCC AR6)

    # ── Economic analysis ─────────────────────────────────────────────────────
    "project_life_yr":  20,         # project lifetime (years)
    "capex_pv_usd_Wp":  0.63,       # APV elevated structure (IRENA 2023 + 40% premium)
    "capex_ad_usd_m3":  800.0,      # AD digester installed cost (Lyberatos 2021)
    "capex_agri_usd_ha":2500.0,     # irrigation + soil preparation (FAO benchmark)
    "bop_fraction":     0.15,       # balance of plant — 15% of PV + AD CAPEX
    "opex_fraction":    0.015,      # annual O&M — 1.5% of total CAPEX

    # ── Digestate N credit (Scenario S4) ─────────────────────────────────────
    "N_content_kg_per_kg_residue": 0.003,   # kg N per kg fresh residue
    "N_price_usd_kg":   1.20,       # nitrogen credit (USD/kg N)

    # ── Biomethane premium (Scenario S5) ─────────────────────────────────────
    "biomethane_premium_factor": 5.5,   # p_biogas = 5.5 × p_elec (EEG 2023 / RED III)
    # NOTE: S5 is labeled "EU policy scenario" — applies directly to Almeria
    # and Freiburg. For Konya and Ouagadougou, use S6 (carbon credit) instead.

    # ── LPG displacement (Scenario S3) ───────────────────────────────────────
    "p_LPG_usd_kWh":    0.18,       # LPG equivalent price (USD/kWh)
}

# =============================================================================
# 4. DIFFERENTIAL EVOLUTION OPTIMIZER SETTINGS
# =============================================================================

DE_SETTINGS = {
    # Design variable bounds — (lower, upper) — order must match VARIABLE_NAMES
    "bounds": [
        (10.0,  40.0),      # beta — panel tilt angle (degrees)
        (4.0,   14.0),      # d_row — row spacing (m)
        (2.0,   5.0),       # H_m — module bottom-edge height (m)
        (2.0,   50.0),      # V_dig — digester volume (m³) — small farm scale
        (20.0,  40.0),      # HRT — hydraulic retention time (days)
        (0.05,  0.50),      # f_PV_heat — PV heating fraction
    ],
    "variable_names": [
        "beta_deg", "d_row_m", "H_m_m",
        "V_dig_m3", "HRT_days", "f_PV_heat"
    ],

    # DE hyperparameters (Storn & Price 1997)
    "population_size":  90,         # NP = 15 × n_variables
    "mutation_F":       (0.5, 1.0), # dithered mutation factor range
    "crossover_Cr":     0.70,       # binomial crossover probability
    "max_generations":  50,
    "seed":             42,         # for reproducibility

    # Feasibility constraints
    "OLR_min":          1.5,        # kg VS/m³·d — minimum organic loading rate
    "OLR_max":          5.0,        # kg VS/m³·d — maximum organic loading rate
    "infeasibility_penalty": -1e6,  # objective value for infeasible solutions

    # Local polishing after DE
    "polish":           True,       # apply L-BFGS-B after DE convergence
    "polish_method":    "L-BFGS-B",

    # NSGA-II settings (water trade-off — Konya only)
    "nsga2_population": 50,
    "nsga2_generations": 80,
}

# =============================================================================
# 5. AD VALORIZATION SCENARIOS (S0–S6)
#    S6 is new — carbon credit scenario for non-EU sites
# =============================================================================

SCENARIOS = {
    "S0": {
        "label":        "Baseline (tomato residues only, grid export)",
        "BMP_mult":     1.00,
        "eta_cap":      1.00,
        "biogas_price": "default",      # uses site p_biogas_usd_kWh
        "digestate_N":  False,
        "carbon_credit":False,
    },
    "S1": {
        "label":        "Co-digestion: residues + cattle manure (2 heads/ha)",
        "BMP_mult":     1.15,
        "eta_cap":      1.00,
        "biogas_price": "default",
        "digestate_N":  False,
        "carbon_credit":False,
    },
    "S2": {
        "label":        "Sub-optimal CH4 capture (85%)",
        "BMP_mult":     1.00,
        "eta_cap":      0.85,
        "biogas_price": "default",
        "digestate_N":  False,
        "carbon_credit":False,
    },
    "S3": {
        "label":        "LPG displacement (on-site use)",
        "BMP_mult":     1.00,
        "eta_cap":      1.00,
        "biogas_price": "LPG",          # p_LPG_usd_kWh from PARAMS
        "digestate_N":  False,
        "carbon_credit":False,
    },
    "S4": {
        "label":        "Digestate N-P-K credit",
        "BMP_mult":     1.00,
        "eta_cap":      1.00,
        "biogas_price": "default",
        "digestate_N":  True,
        "carbon_credit":False,
    },
    "S5": {
        "label":        "Biomethane grid injection (EU EEG/RED III premium)",
        "BMP_mult":     1.00,
        "eta_cap":      1.00,
        "biogas_price": "biomethane",   # 5.5 × p_elec
        "digestate_N":  False,
        "carbon_credit":False,
        "note":         "EU policy scenario — directly applicable to Almeria and Freiburg",
    },
    "S6": {
        "label":        "Carbon credit (voluntary market)",
        "BMP_mult":     1.00,
        "eta_cap":      1.00,
        "biogas_price": "default",
        "digestate_N":  False,
        "carbon_credit":True,           # uses site carbon_price_usd_tCO2
        "note":         "Non-EU scenario — applicable to Konya and Ouagadougou",
    },
}

# =============================================================================
# 6. CSV OUTPUT COLUMN DEFINITIONS
#    Documents every column that will appear in results_summary.csv
#    so the structure is decided before writing any simulation code.
# =============================================================================

CSV_COLUMNS = [
    # Identifiers
    "site", "scenario",

    # Optimal design variables (from optimizer)
    "beta_deg", "d_row_m", "H_m_m", "V_dig_m3", "HRT_days", "f_PV_heat",
    "GCR_pct", "OLR_kgVS_m3d",

    # E1 — Energy outputs (annual, per hectare)
    "PV_total_MWh_ha", "PV_sold_MWh_ha", "PV_self_MWh_ha",
    "biogas_total_MWh_ha", "heat_demand_MWh_ha",
    "biogas_for_heat_MWh_ha", "ESR",
    "BMP_avg_NmL_gVS",

    # eLER decomposition
    "LER_crop", "LER_PV_defA", "LER_PV_defB", "LER_biogas",
    "eLER_defA", "eLER_defB", "LER_2C",

    # Water savings — CORRECTED: season denominator
    "W_saved_mm_season",         # mm saved over growing season
    "ET0_season_mm",             # growing-season ET0 (denominator)
    "ET0_annual_mm",             # annual ET0 (for reference only)
    "LER_water_season_pct",      # W_saved / ET0_season × 100  ← paper value
    "LER_water_annual_pct",      # W_saved / ET0_annual × 100  ← old (wrong) value

    # E2 — Economic outputs
    "CAPEX_PV_kUSD", "CAPEX_AD_kUSD", "CAPEX_agri_kUSD",
    "CAPEX_BOP_kUSD", "CAPEX_total_kUSD",
    "gross_revenue_kUSD_yr", "opex_kUSD_yr", "net_CF_kUSD_yr",
    "NPV_kUSD", "IRR_pct", "payback_yr",

    # E3 — Energo-economic outputs
    "LCOE_PV_USD_kWh",           # levelised cost of PV electricity
    "cost_biogas_USD_MWh",       # levelised cost of biogas production

    # E4 — Environmental outputs
    "CO2_avoided_PV_tCO2_ha",    # from grid electricity displacement
    "CO2_avoided_biogas_tCO2_ha",# from natural gas displacement
    "CO2_avoided_total_tCO2_ha",
    "carbon_value_USD_ha",       # enviro-economic: CO2_total × carbon_price

    # Validation flags
    "V1_MAPE_pct",               # crop yield model MAPE (filled in 06_validation)
    "V2_MAPE_pct",               # BMP model MAPE
]

# =============================================================================
# 7. FIGURE SETTINGS
#    Consistent style across all 8 figures
# =============================================================================

FIG_STYLE = {
    "dpi":          300,            # publication quality
    "format":       "png",
    "figsize_single": (10, 6),      # single-panel figure
    "figsize_double": (14, 6),      # two-panel figure
    "figsize_quad":   (14, 10),     # four-panel figure (one per site)

    # Color palette — one color per site, consistent across all figures
    "site_colors": {
        "Konya":        "#1f77b4",  # blue
        "Almeria":      "#d62728",  # red
        "Ouagadougou":  "#2ca02c",  # green
        "Freiburg":     "#8c510a",  # brown
    },

    # eLER component colors — consistent with seminar presentation
    "eLER_colors": {
        "LER_crop":   "#4c87c2",    # blue
        "LER_PV":     "#f5a623",    # orange
        "LER_biogas": "#2e9b5a",    # green
    },

    "benchmark_LER": 1.94,          # Riaz et al. 2022 — red dashed line on eLER figures
}

# =============================================================================
# 8. QUICK SELF-CHECK
#    Run this file directly to verify all paths and parameters are consistent.
#    python 00_config.py
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("APV-AD Project — Configuration check")
    print("=" * 60)

    # Check data files
    print("\n[1] PVGIS data files:")
    all_found = True
    for site, path in PVGIS_FILES.items():
        status = "✓ found" if path.exists() else "✗ NOT FOUND"
        if not path.exists():
            all_found = False
        print(f"    {site:15s} → {path.name:35s} {status}")

    # Check output directories
    print("\n[2] Output directories:")
    for d in [CSV_DIR, FIG_DIR, PAPER_DIR]:
        print(f"    {'✓' if d.exists() else '✗'} {d}")

    # Check DE bounds vs variable names
    n_bounds = len(DE_SETTINGS["bounds"])
    n_names  = len(DE_SETTINGS["variable_names"])
    de_ok    = n_bounds == n_names == 6
    print(f"\n[3] DE settings: {n_bounds} bounds, {n_names} variable names"
          f" — {'✓ consistent' if de_ok else '✗ MISMATCH'}")

    # Check scenarios
    print(f"\n[4] Scenarios defined: {list(SCENARIOS.keys())}")

    # Check CSV columns
    print(f"\n[5] CSV output columns: {len(CSV_COLUMNS)} columns defined")

    # Check sites
    print(f"\n[6] Sites configured: {list(SITES.keys())}")
    for site, cfg in SITES.items():
        vs_residues = (cfg["fruit_yield_t_ha"] * 1000
                    * cfg["R_res"] * cfg["DM"] * cfg["VS_fraction"])
        vs_manure   = (cfg["manure_input_t_ha"] * 1000
                       * cfg["manure_DM"] * cfg["manure_VS_fraction"])
        vs_total    = vs_residues + vs_manure
        print(f"    {site:15s} → VS_residues={vs_residues:.0f} kg  VS_manure={vs_manure:.0f} kg  VS_total={vs_total:.0f} kg/ha/yr"
              f"  |  discount = {cfg['discount_rate']*100:.0f}%"
              f"  |  ef_grid = {cfg['grid_ef_tCO2_MWh']} tCO2/MWh")

    print("\n" + "=" * 60)
    if all_found and de_ok:
        print("All checks passed. Ready to run 01_load_data.ipynb")
    else:
        print("WARNING: Fix the issues above before proceeding.")
    print("=" * 60)
