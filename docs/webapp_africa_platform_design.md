# Africa APV+AD Platform — Design Brainstorm & Integration Plan

Status: brainstorm / pre-implementation design notes. Not yet built.
Companion to the existing validation docs in this folder (`pvlib_validation.md`,
`ad_bmp_sensitivity.md`, `eLER_weighting_sensitivity.md`, `et0_ra_fix_notes.md`,
`poa_azimuth_fix_notes.md`).

## 1. Concept

A public web platform with two layers:

1. **Continental potential maps** — wind and solar resource choropleths across
   Africa (Global Wind Atlas / PVGIS-Global Solar Atlas as data sources), the
   "explore" entry point.
2. **Point-and-analyze APV+AD presizing tool** — click any location, pick a
   crop and co-substrate, and get an eLER decomposition + 4E dashboard
   (NPV/IRR/LCOE/CO2) computed by the *actual validated model from this repo*,
   not a simplified web reimplementation.

This is not a new model — it is `notebooks/simulation_functions.py` +
`notebooks/config_00.py` (the exact six-model pipeline behind the manuscript)
wrapped in a web API and generalized from "4 fixed sites" to "any point in
Africa." The manuscript's own four sites (Konya, Almería, Ouagadougou,
Freiburg) become built-in, pre-validated reference cases inside the tool —
Ouagadougou in particular, since it's the only African site already studied
and gives a direct trust anchor ("run this tool on Ouagadougou, get Table 8/9/12
back exactly").

Target audience: researchers/students (per discussion) — so the tool should
expose methodological transparency (Definition A/B toggle, weighting scheme,
confidence tiers, validation stats) rather than hide it behind a simplified UI.

## 2. What already exists vs. what needs building

### Already exists (reuse, don't rebuild)

- `notebooks/simulation_functions.py`: `compute_solar_angles`, `compute_shading`,
  `compute_POA`, `compute_PV`, `compute_crop_yield`, `compute_ET_reduction`,
  `compute_AD`, `eLER_weighted`, `eLER_objective`, `modules_per_ha`,
  `reference_pv_kWh` — the full validated pipeline, pytest-covered.
- `notebooks/config_00.py`: `SITES`, `PARAMS`, `DE_SETTINGS`, `SCENARIOS`,
  `CSV_COLUMNS` — single source of truth for constants.
- `notebooks/03_optimizer.ipynb`: DE + L-BFGS-B optimization workflow
  (scipy.optimize.differential_evolution), matching Table 5/8 exactly.
- `notebooks/07_pareto_analysis.py`: NSGA-II bi-objective (eLER vs. LER_water)
  script, already runnable standalone with `--site`/`--pop`/`--gen` flags.
- `eLER_weighted(w_crop, w_PV, w_bio)`: weighting-scheme sensitivity, already
  implemented and backward-compatible (default reproduces unweighted eLER).
- `compute_AD(..., AD_overrides=...)`: Monte Carlo sensitivity hook for
  `BMP_ref`, capture efficiency, VS fractions — already implemented.
- `tests/test_simulation_functions.py`, `tests/test_notebook_smoke.py`: pytest
  suite parametrized across all 4 sites, plus a real-PVGIS-data smoke test.

### Needs building

- **Site-generalization refactor**: `compute_AD(df, site: str, ...)` and
  `eLER_objective(x, df, site: str, ...)` currently resolve `site` via
  `SITES[site]` lookup. Change the type to `site: str | dict` and resolve with
  `cfg = SITES[site] if isinstance(site, str) else site`. Minimal,
  backward-compatible (existing string-keyed calls in the notebooks/tests are
  unaffected), and is the single change that unlocks "any African point" —
  the API layer builds a one-off dict per request instead of requiring a
  registered site name.
- **Data ingestion service**: replicate `01_load_data.ipynb`'s TMY→DataFrame
  pipeline (`GHI, DNI, DHI, T2m, WS, in_season` columns) for an arbitrary
  `(lat, lon)` via the live PVGIS API, instead of reading the four checked-in
  CSVs. Cache results (PVGIS TMY doesn't change often).
- **Crop/substrate preset tables**: see Section 3 — same field names as
  `SITES[site]`'s agronomic block, extended with a `confidence_tier`.
- **Growing-season determination**: `SITES[site]["planting_date"/"harvest_date"]`
  are hand-set per site today. For arbitrary locations, default to an
  FAO/GIEWS crop-calendar lookup by (country, crop), with a climate-inferred
  fallback and user override — never silently invented.
- **Web API** (FastAPI) wrapping the above — see Section 4.
- **Hemisphere check**: `compute_solar_angles`/`compute_POA` assume a
  due-south-facing fixed-tilt collector (correct for all 4 Northern-Hemisphere
  study sites). Africa spans the equator substantially — verify/adjust for
  Southern Hemisphere points (north-facing optimal orientation) before trusting
  results south of ~5°N. This is a code-correctness check, not a data-sourcing
  gap, and none of the 4 original sites would have exposed it.
- **`min_beta_deg` generalization**: each `SITES[site]` entry hardcodes a
  latitude-based minimum tilt (20°, 20°, 12°, 22°) rather than deriving it from
  a formula. Needs a general rule (paper's Finding 2: optimal tilt ≈
  0.4–0.5 × latitude) before it can be computed for an arbitrary point instead
  of looked up.

## 3. Preset schema (maps directly onto existing `SITES[site]` fields)

```python
# Crop presets — extends the agronomic block already in SITES[site]
crop_presets = {
    "tomato": {
        "PAR_sat": 174.0, "fruit_yield_t_ha": 60.0, "R_res": 0.80,
        "DM": 0.06, "VS_fraction": 0.85,
        "confidence_tier": "validated",   # V1 in 06_validation.ipynb
        "source": "Mohammedi et al. 2023",
    },
    "lettuce":  {"PAR_sat": 76.0,  "confidence_tier": "literature_estimated", ...},
    "spinach":  {"PAR_sat": 80.0,  "confidence_tier": "literature_estimated", ...},
    "maize":    {"PAR_sat": 400.0, "confidence_tier": "literature_estimated", ...},
    # residue_ratio / DM / yield for non-tomato crops: NOT filled in yet —
    # placeholder gap, needs real per-crop sourcing before use (see Section 6).
}

# Substrate presets — extends manure_* fields already in SITES[site]
substrate_presets = {
    "cattle_manure": {
        "DM": 0.20, "VS_fraction": 0.80, "BMP_ref": 300.0,
        "confidence_tier": "validated",   # Pilarski 2025, used throughout paper
    },
    "poultry_litter": {"DM": 0.28, "VS_fraction": 0.72, "confidence_tier": "literature_estimated"},
    "goat_manure":    {"DM": 0.32, "VS_fraction": 0.78, "confidence_tier": "literature_estimated"},
    # BMP_ref left unset for non-cattle substrates pending real citations —
    # do not guess a methane yield number, it feeds directly into 4E economics.
}
```

`confidence_tier` values: `validated` (has its own field/lab dataset like the
paper's V1/V2), `literature_estimated` (single-source value, no
temperature/shading response curve), `user_custom` (no citation). Anything
below `validated` should auto-trigger the `AD_overrides`-based Monte Carlo
that already exists in `compute_AD`, and surface the resulting range in the UI
rather than a bare point estimate.

## 4. API sketch

```
GET  /api/climate?lat=&lon=
  → fetch/cache PVGIS TMY, build the same DataFrame schema as 01_load_data,
    return summary (annual GHI, mean T) for map tooltips

POST /api/evaluate        (fast tier — single eLER_objective call, no DE)
  body: { lat, lon, crop, co_substrates: [...], scenario,
          design_vars?, econ_overrides?, weights?: {w_crop,w_PV,w_bio} }
  → builds a runtime site-dict (per the str|dict refactor), calls
    eLER_objective(x, df, site=site_dict, return_full=True, ...)

POST /api/optimize        (full tier — background job)
  body: same as /evaluate minus design_vars
  → runs the same DE + L-BFGS-B as 03_optimizer.ipynb, streams
    generation/best-eLER progress (reproduces Fig. 3 convergence chart live)

POST /api/pareto           (background job)
  → wraps 07_pareto_analysis.py's NSGA-II search for the given site-dict

GET  /api/sites/reference  → the 4 published sites' exact Table 8/9/12 values,
                              pulled from outputs/csv/results_summary.csv,
                              for the trust-anchor comparison
```

Correctness gate before anything else ships: `/api/evaluate` called on
Ouagadougou's exact `(lat=12.37, lon=-1.53)` with tomato + cattle manure + S0
must reproduce `results_summary.csv`'s Ouagadougou row exactly. If it doesn't,
nothing downstream is trustworthy.

## 5. Suggested build order

1. `site: str | dict` refactor in `simulation_functions.py` (small, tested
   against existing 4-site pytest suite to confirm no regression).
2. Data ingestion service (PVGIS API → DataFrame), validated by re-deriving
   the four existing `data/*_climate.csv` files from live PVGIS and diffing.
3. Thin FastAPI wrapper exposing `/api/evaluate` only, tomato + cattle manure
   only, hardcoded scenario S0 — prove Ouagadougou reproduction (Section 4's
   correctness gate) before anything else.
4. `/api/optimize` (DE job + live convergence stream).
5. Map layer (wind/solar choropleth) as the front door into `/api/evaluate`.
6. Multi-crop / multi-substrate presets — **only after** each new preset has
   either its own validation dataset or an explicit, visible Monte Carlo range
   (see Section 6). Do not ship placeholder PAR_sat/BMP_ref values as if they
   were validated tomato/cattle-manure-grade numbers.
7. `/compare` page, weighting-scheme toggle UI, confidence-tier badges.

## 6. Generalization limits — gates before expanding beyond tomato/cattle-manure/4-sites

- **New crop**: needs its own PAR-response validation (like V1's shading-trial
  MAPE check), not just a PAR_sat literature value — the yield model's hard
  PAR-saturation threshold is a specific functional form validated for tomato
  under Mediterranean greenhouse shading (Fshad ≤ 0.25), and C4 crops (maize,
  sorghum, sugarcane) may need a different functional form entirely, not just
  a different threshold.
- **New substrate**: needs its own multi-temperature BMP dataset (like V2's
  Feng et al. curve for tomato residues) before trusting the shared Arrhenius
  θ=0.069°C⁻¹ — that coefficient encodes tomato-residue microbial kinetics
  specifically. No inhibition model exists yet (only the OLR feasibility
  bound); arbitrary substrate mixing needs at least a C:N sanity check before
  the tool allows it.
- **New location**: boundary-convergence (row spacing/height at physical
  minimum) is confirmed only across the paper's GHI 1150–2050 kWh/m² and
  11.4–28.5°C envelope. Hargreaves-Samani ET0 ignores humidity and was
  validated against Penman-Monteith for semi-arid/temperate climates
  specifically — humid tropical zones (Congo Basin, coastal West Africa) are
  exactly where that omission is expected to matter most and should be
  spot-checked, not assumed fine.
- **Economics**: `SITES[site]`'s CAPEX constants (`capex_pv_usd_Wp`,
  `capex_ad_usd_m3`) are calibrated to European/Turkish supply chains. A
  landlocked-country CAPEX multiplier (IRENA Africa cost reports) should be
  applied rather than reusing one global constant continent-wide.

General rule: nothing graduates out of `literature_estimated`/`user_custom`
tier into a confident headline number until it has either a real validation
dataset or a visible Monte Carlo range — a polished dashboard on top of an
unvalidated assumption is worse than an obviously incomplete one.
