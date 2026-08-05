# APV + AD Code Analysis

## 📋 Project Overview

**Integrated Agrivoltaic (APV) – Anaerobic Digestion (AD) Optimization Model**
- **Version**: 3.0 (Production Release)
- **Status**: Ready for journal submission (Renewable & Sustainable Energy Reviews / Applied Energy)
- **Study Sites**: 4 Köppen climate zones across 3 continents

### Study Sites Configuration

| Site | Location | Climate | Coordinates | Latitude | Longitude |
|------|----------|---------|-------------|----------|-----------|
| **Konya** | Turkey | BSk (Cold semi-arid) | 37.87°N 32.49°E | 37.87 | 32.49 |
| **Almería** | Spain | BSh (Hot semi-arid) | 36.83°N 2.46°W | 36.83 | -2.46 |
| **Ouagadougou** | Burkina Faso | BSh (Tropical semi-arid) | 12.37°N 1.52°W | 12.37 | -1.52 |
| **Freiburg** | Germany | Cfb (Oceanic temperate) | 47.99°N 7.85°E | 47.99 | 7.85 |

**Crop**: Tomato (Solanum lycopersicum)
- Semi-arid sites: 60 t/ha irrigated yield
- Freiburg: 30 t/ha open-field yield (DLG, 2023)

---

## 🔧 Key Files and Their Purposes

### 1. **APV_AD_v3.py** (88.5 KB) — Main Simulation Engine
**Core Functionality:**
- Integrated bifacial agrivoltaic PV system modeling
- Anaerobic digestion energy valorization
- Four-site comparative analysis
- Multi-objective optimization with Pareto frontiers

**Key Classes:**

#### SolarGeometry (Lines 255–331)
- **Purpose**: Hour-by-hour solar position calculations
- **Methods**:
  - `solar_declination()`: Spencer (1971) declination formula
  - `equation_of_time()`: True solar time correction
  - `solar_altitude()`, `solar_azimuth()`: Position angles
  - `shading_factor()`: Ground shading (Zainali et al., 2023)
- **Inputs**: Latitude, longitude, timezone offset
- **Precision**: Sub-degree accuracy for row orientation, tilt angle optimization

#### AgrivoltaicSystem (Lines 337–456)
- **Purpose**: Bifacial PV system + crop interaction modeling
- **Physical Model**:
  - **Module**: Jinko Tiger Neo N-type bifacial, 550 Wp, bifaciality = 0.70
  - **Irradiance**: Liu & Jordan (1960) plane-of-array method
  - **Cell Temperature**: Faiman (2008) model (U0=25 W/m²·K, U1=6.84 W·s/m³·K)
  - **Crop PAR Yield**: **FIX 1 (Critical)** — Asymmetric PAR saturation
    - Denominator: Raw PAR_open integral (NO saturation clip)
    - Numerator: PAR_AV clipped at 174 W/m² (PAR_sat)
    - Restores LER_crop to physically correct range 0.55–0.65
    - Fixes LER_crop = 1.000 bug from earlier versions
- **Evapotranspiration**: FAO-56 hourly Penman-Monteith (site-specific γ)
- **Parameters Varied in R1 (Monte-Carlo)**:
  - `alpha_shade`: Shading ET reduction coefficient

#### AnaerobicDigestionSystem (Lines 462–500+)
- **Purpose**: Biomass-to-energy conversion via mesophilic AD
- **BMP Model**: Temperature-corrected Arrhenius (Pilarski & Pilarska, 2025)
- **Parameters R1 (varied)**:
  - `BMP_ref`: Reference BMP at 37°C (NmL CH4/g VS)
  - `U_eff`: Digester wall U-value (W/m²·K)
  - `fruit_yield_t_ha`: Site-specific (FIX 4)

---

### 2. **apv_ad_pv_validation.py** (392 lines) — PVGIS Reference Validation
**Purpose**: Validate simulation PV outputs against PVGIS-SARAH3 calculator

**Key Features:**
- Calls PVGIS PVcalc API for 4 sites at optimal tilts (Table 9)
- Compares annual yield (kWh/kWp) with simulation
- Accounts for bifacial rear-gain bias (~8.4%)
- Generates:
  - Console validation table
  - LaTeX Table V3 for paper
  - JSON cache for offline runs

**Validation Metrics:**
- **Raw MAPE**: Includes bifacial rear gain
- **Adjusted MAPE**: Bifacial gain removed (~8.4%)
- **Expected Range**: Within PVGIS-SARAH3 ±5% monthly GHI uncertainty

**PVGIS Parameters:**
- System: 1 kWp c-Si, loss=12%, south-facing
- Database: PVGIS-SARAH3
- Technology: crystSi (monocrystalline)

---

## 🔨 Version 3.0 Fixes and Improvements

### CRITICAL FIXES

**FIX 1** — **crop_yield_ratio() asymmetric PAR saturation** ✅
- **Issue**: Previous versions clipped PAR_open denominator, forcing LER_crop → 1.000
- **Root Cause**: Symmetric saturation violated physical reality; semi-arid open fields exceed PAR saturation at midday
- **Solution**: 
  - Denominator = raw PAR_open integral (no clip)
  - Numerator = PAR_AV clipped at PAR_sat only
- **Impact**: LER_crop now 0.55–0.65 (physically correct) vs. spurious 1.000

**FIX 2** — LER_PV reference cache with tilt invalidation
- Tracks tilt used in `_ref_pv_tilt_cached`
- Rebuilds cache when tilt changes (optimization sweeps)
- Prevents stale reference PV values across tilt scenarios

**FIX 3** — Ouagadougou cross-year growing season
- Corrected: Jul 1 – Oct 30 (was Oct 1 – Oct 30)
- Reflects wet monsoon cycle properly
- OR-logic guard ensures robustness against future cross-year sites

**FIX 4** — Site-specific economics & crop yield
- p_elec, p_biogas, p_tomato, fruit_yield_t_ha stored per-site
- No global ECONOMICS mutation
- Enables realistic 30 t/ha (Freiburg) vs. 60 t/ha (semi-arid)

**FIX 5** — AD Scenario S5 (biomethane grid injection)
- Valued at **5.5 × p_elec** (Nik Zad et al., 2025; EEG 2023)
- Reflects German grid-injection premium for renewables-based biomethane
- Single SCENARIOS entry; no structural change to energy balance

**FIX 6** — Piecewise BMP temperature model (validation-only)
- **T ≤ 37°C**: Arrhenius model (Pilarski & Pilarska, 2025)
- **T > 37°C**: Linear inhibition −5 NmL/gVS/°C (Feng et al., 2013)
- T_eff clipped at 37°C operationally → inhibition branch used only in `run_validation()`

**FIX 7** — Monte-Carlo rewritten as structural invariance + economic sensitivity
- **Replaces**: CV=0.0% spike histograms
- **New Approach**:
  - (a) One-way tornado on design variables (genuine eLER sensitivity)
  - (b) Net income P5/P95 distribution (genuine economic Monte-Carlo output)

**FIX 8** — Freiburg site added
- Location: 47.99°N, 7.85°E, Cfb
- Yield: 30 t/ha (open-field German baseline)
- Season: May 1 – Oct 15
- γ = 0.0673 kPa/°C at 278 m elevation

**FIX 9** — Pareto base class safe under missing pymoo
- Runtime-resolved base class avoids NameError if pymoo not installed
- Graceful degradation with informative message

**FIX 10** — Pareto site selection (season-length gate)
- Pareto runs only if growing season > 150 days
- Konya: 189d ✅, Freiburg: ~167d ✅
- Almería: 122d ❌, Ouagadougou: 122d ❌
- Excludes sites where LER_water variance below numerical resolution

### PERFORMANCE IMPROVEMENTS

**PERF 1** — Vectorized reference PV computation
- Replaced row-by-row `iterrows()` loop with full numpy vectorization
- **Speedup**: ~400× faster

**PERF 2** — Monte-Carlo geometry pre-computation
- Reuse pre-computed `solar_altitude_rad`, `solar_azimuth_rad` arrays
- Only vary physics parameters (BMP, U_eff, α_shade) and prices
- **Speedup**: ~30× for N=2000 samples

---

## 📊 Configuration and Constants

### RUN_CONFIG Dictionary
```python
RUN_CONFIG = {
    "sites": ["Konya", "Almeria", "Ouagadougou", "Freiburg"],
    "de_maxiter": 50,
    "de_popsize": 12,  # → pop = 72
    "de_seed": 42,
    "run_monte_carlo": True,
    "mc_n_samples": 2000,
    "mc_seed": 99,
    "run_mo_constraint_sweep": True,
    "run_mo_pareto": True,
    "mo_d_row_values": [3, 4, 5, 6],  # row spacing (m)
    "pareto_min_season_days": 150,
    "run_ad_scenarios": True,
    "run_validation": True,
    "save_figures": True,
    "figure_dpi": 300,
    "output_dir": "apv_ad_results",
}
```

### CAPEX & Economics
```python
CAPEX = {
    "pv_capex_usd_module": 280.0,        # Jinko Tiger Neo 550 Wp (2024)
    "ad_capex_usd_m3": 450.0,            # farm-scale mesophilic digester
    "struct_capex_frac": 0.40,           # mounting structure
    "opex_frac_capex": 0.025,            # annual O&M
    "p_N_fertilizer_usd_kg": 1.20,       # N-equivalent for digestate (S4)
    "lpg_price_usd_kWh": 0.18,           # LPG parity for local use (S3)
}
BIOMETHANE_PREMIUM = 5.5  # × p_elec for scenario S5
```

### AD Valorization Scenarios
| Scenario | Label | BMP | Capture | Local | Biomethane | Digestate | Key Parameter |
|----------|-------|-----|---------|-------|-----------|-----------|----------------|
| **S0** | Baseline (mono-dig, grid export) | 1.0× | 1.0× | ❌ | ❌ | ❌ | Baseline reference |
| **S1** | Co-digestion (30% manure, +15% BMP) | 1.15× | 1.0× | ❌ | ❌ | ❌ | BMP_mult=1.15 |
| **S2** | Sub-optimal (85% CH4 capture) | 1.0× | 0.85× | ❌ | ❌ | ❌ | eta_cap=0.85 |
| **S3** | Local LPG displacement | 1.0× | 1.0× | ✅ | ❌ | ❌ | on-site use |
| **S4** | Digestate valorization (N-P-K) | 1.0× | 1.0× | ❌ | ❌ | ✅ | fertilizer credit |
| **S5** | Biomethane grid injection | 1.0× | 1.0× | ❌ | ✅ | ❌ | 5.5× p_elec |

---

## 🧮 Physical Models & Equations

### 1. Solar Geometry (Spencer, 1971)

**Declination (δ):**
$$\delta = 23.45° \sin\left(\frac{360°}{365}(n - 81)\right)$$

**Hour Angle (ω):** 15° per solar hour
$$\omega = 15° \times (t_{solar} - 12:00)$$

**Solar Altitude (α_s):**
$$\sin(\alpha_s) = \sin(\phi) \sin(\delta) + \cos(\phi) \cos(\delta) \cos(\omega)$$

### 2. Plane-of-Array Irradiance (Liu & Jordan, 1960)

$$G_{POA,front} = G_{beam} \cos(\theta) + G_{diff} \frac{1+\cos(\beta)}{2} + G_{albedo} \frac{1-\cos(\beta)}{2}$$

**Bifacial Effective (with φ=0.70 bifaciality):**
$$G_{eff} = G_{front} + φ \times G_{rear}(1 - \eta_{loss})$$

### 3. Cell Temperature (Faiman, 2008)

$$T_{cell} = T_{amb} + \frac{G_{POA}}{U_0 + U_1 \cdot \max(0.5, v_{wind})}$$

Where:
- U₀ = 25.0 W/m²·K (constant thermal loss)
- U₁ = 6.84 W·s/m³·K (wind-dependent)

### 4. PV Power Output (Eq. 10)

$$P = n_{modules} \times P_{STC} \times \frac{G_{eff}}{1000} \times [1 - \eta_{loss}] \times [1 + \gamma_{Pmax}(T_{cell} - 25°C)]$$

Where:
- γ_Pmax = -0.0035 /°C (Jinko temp coefficient)
- η_loss = 12% (wiring, inverter, mismatch)

### 5. Crop Yield Ratio — FIX 1 (Asymmetric PAR Saturation)

$$LER_{crop} = \frac{\int \min(PAR_{AV}, PAR_{sat}) \, dt}{\int PAR_{open} \, dt}$$

**Critical Correction:**
- Numerator: PAR_AV clipped at PAR_sat = 174 W/m²
- **Denominator: PAR_open NOT clipped (raw integral)**

**PAR Relations:**
- PAR_open = GHI × 0.48 (McCree, 1971)
- PAR_AV = PAR_open × (1 − F_shad)
- PAR_sat = 174 W/m² (800 µmol/m²/s ÷ 4.6 µmol/W)

### 6. FAO-56 Hourly Penman-Monteith ET₀ (Eq. 12)

$$ET_0 = \frac{0.408 \Delta R_n + \gamma \frac{37}{T+273} u_s (e_s - e_a)}{\Delta + \gamma(1 + 0.24 u_s)}$$

Where:
- Δ = slope of saturation vapor pressure
- γ = psychrometric constant (site-specific, 0.0665–0.0673 kPa/°C)
- R_n ≈ 0.77 × GHI (net radiation)

### 7. APV ET Reduction (Eq. 14)

$$ET_{APV} = ET_0 \times [1 - \alpha_{shade} \times F_{shad}]$$

- α_shade ∈ [0.20, 0.60] (varied in Monte-Carlo R1)
- F_shad ∈ [0, 1] from solar geometry

### 8. Anaerobic Digestion BMP (Temperature-Corrected)

**At T ≤ 37°C (Arrhenius, Pilarski & Pilarska 2025):**
$$BMP(T) = BMP_{ref} \times \exp\left[\alpha_{BMP} (T - 37)\right]$$

**At T > 37°C (Linear Inhibition, Feng et al. 2013):**
$$BMP(T) = BMP_{ref} - 5 \times (T - 37)$$

- BMP_ref ∈ [250, 350] NmL CH₄/g VS (R1 parameter)
- Reference: 37°C (optimal mesophilic)
- Validation uses piecewise model; operations use T_eff clipped at 37°C

### 9. Biomass-to-Biogas Energy

$$E_{biogas} = \text{Residue Mass (t)} \times \text{BMP}(T) \times \text{CH4 LHV}$$

- Residue from: crop yield × specific waste factor
- LHV_CH₄ = 35.8 MJ/m³ STP
- η_capture ∈ {0.85 (S2), 1.00 (others)}

---

## 📈 Analysis Workflows

### R1 — Monte-Carlo Sensitivity
**Varied Parameters (N=2000 samples):**
1. BMP_ref: [250, 350] NmL CH₄/g VS
2. U_eff: [0.5, 1.2] W/m²·K
3. α_shade: [0.20, 0.60] (ET reduction)

**Outputs:**
- Tornado diagram (eLER sensitivity)
- Net income P5/P95 distribution
- System robustness assessment

### R2 — Constraint Sweep (d_row ∈ {3, 4, 5, 6} m)
**Varied:** Row spacing d_row
**Fixed:** All other parameters at baseline
**Outputs:**
- LER vs. row spacing curves
- Optimal trade-off zones
- GCR (ground cover ratio) ranges

### R3-A — Single-Objective Optimization (DE)
**Algorithm:** Differential Evolution (scipy.optimize)
- max_iter = 50
- pop = 12 × n_var = 72
- seed = 42

**Objective:** Maximize eLER (economic LER)
- Accounts for both biomass value & PV revenue

**Variables:** Row spacing, tilt angle, AD design parameters

### R3-B — Multi-Objective Pareto (NSGA-II)
**Algorithm:** NSGA-II (pymoo)
**Sites:** Konya, Freiburg only (growing season > 150d)
**Objectives:**
1. Maximize LER_crop (crop yield preservation)
2. Maximize LER_water (water sustainability)
3. Maximize NPV or eLER (economic viability)

**Output:** Pareto frontier visualizations

---

## 📁 Supporting Data & Outputs

### Input Files
- `konya_climate.csv` — PVGIS-SARAH3 TMY (hourly)
- `almeria_climate.csv`
- `ouagadougou_climate.csv`
- `freiburg_climate.csv`

**Format per file:** Year, Month, Day, Hour, GHI, DNI, DHI, T_amb, RH, WS

### Output Directory: `apv_ad_results/`
- `*.png` — Figures (300 dpi)
- `*.pdf` — Publication-ready plots
- `summary.csv` — Site & scenario results table
- `table_v3_pv_validation.tex` — LaTeX Table V3
- `discussion_v3_paragraph.tex` — Discussion text for paper

---

## 🐍 Dependencies

### Core
- **numpy** ≥ 1.19.0 — Vectorized computations
- **pandas** ≥ 1.0.0 — Time series climate data
- **scipy** ≥ 1.5.0 — Optimization, statistics
- **matplotlib** ≥ 3.1.0 — Figures & publication plots

### Optional
- **pymoo** ≥ 0.6.0 — NSGA-II Pareto optimization (R3-B)
  - If missing: R3-B disabled, informative warning printed
- **requests** (for apv_ad_pv_validation.py only)
  - PVGIS API calls
  - Offline mode uses cached JSON if requests unavailable

### Installation
```bash
pip install numpy pandas scipy matplotlib pymoo requests
```

---

## 🚀 Usage Quick Start

### 1. Prepare Climate Data
Place PVGIS-SARAH3 TMY CSV files in working directory or update `SITES` paths.

### 2. Configure RUN_CONFIG
```python
RUN_CONFIG = {
    "sites": ["Konya", "Almeria", "Ouagadougou", "Freiburg"],
    "run_monte_carlo": True,
    "run_mo_pareto": True,
    "output_dir": "apv_ad_results",
}
```

### 3. Run Main Simulation
```bash
python APV_AD_v3.py
```

### 4. Validate PV Outputs
```bash
python apv_ad_pv_validation.py              # with internet
python apv_ad_pv_validation.py --offline     # use cached PVGIS
```

### 5. Inspect Results
- Figures in `apv_ad_results/`
- LaTeX tables for paper
- CSV summary: site yields, LER metrics, economics

---

## ✅ Validation Checkpoints

1. **PV Model**: Validated vs. PVGIS-SARAH3 (within ±5% GHI uncertainty)
2. **Crop PAR Saturation**: LER_crop ∈ [0.55, 0.65] ✅
3. **ET Model**: Against FAO-56 reference stations
4. **BMP Model**: Temperature-corrected Arrhenius validated @ 37°C
5. **Economic Parameters**: Cross-checked vs. recent market data (2023–2025)

---

## 📝 Key Paper References

- **Zainali et al. (2023)** — Shading factor model, Appl. Energy 339, 120981
- **Faiman (2008)** — Cell temperature model, Sol. Energy Mater. Cell 92
- **Pilarski & Pilarska (2025)** — BMP temperature correction
- **Feng et al. (2013)** — AD inhibition at high temperature
- **Huld et al. (2012)** — PVGIS-SARAH3 calibration & uncertainty
- **Urraca et al. (2017)** — PVGIS GHI validation
- **Nik Zad et al. (2025)** — Biomethane grid injection premium (EEG 2023)

---

## 🎯 Summary

The APV + AD model (v3.0) is a comprehensive, production-ready system for:
- **Bifacial agrivoltaic PV design** across diverse climates
- **Anaerobic digestion** of crop residues (6 valorization scenarios)
- **Integrated energy & food security analysis**
- **Multi-objective optimization** for trade-offs between crop production, water, and renewable energy

**Key Innovations:**
1. **FIX 1** (asymmetric PAR saturation) — Physically correct crop yield modeling
2. **FIX 4** (site-specific parameters) — Realistic economic & agronomic diversity
3. **FIX 5** (biomethane premium) — Reflects 2025 German renewable-gas market
4. **PERF 1-2** (vectorization & pre-computation) — ~400× faster simulations
5. **Robust validation** against PVGIS reference & peer-reviewed models

**Status**: Ready for peer review & publication.
