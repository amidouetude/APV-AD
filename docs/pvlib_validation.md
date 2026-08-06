# pvlib Validation

This note documents the external validation of the plane-of-array (POA)
irradiance model against `pvlib` (v0.15.2), performed in response to peer
review requesting benchmark comparison against an independent,
widely-used solar modeling library.

## Method

For each site, at its optimized tilt angle (Table 8) and using the site's
actual PVGIS-SARAH3 TMY data (`DNI`, `GHI`, `DHI`):

1. Computed solar position via `pvlib.solarposition.get_solarposition()`
   (NREL SPA algorithm) rather than this study's own Spencer/declination
   formulas, to get an independent ground truth for solar geometry.
2. Computed POA irradiance via `pvlib.irradiance.get_total_irradiance()`
   under two transposition models: `isotropic` (Liu & Jordan, the same
   sky model this study uses — Eq. 7) and `perez` (anisotropic, accounts
   for circumsolar and horizon brightening).
3. Compared annual front-surface POA (`G_front`, Eq. 7) from this study's
   own `compute_POA()` against both `pvlib` outputs.

Ground albedo (0.25) and surface azimuth (180°, due south) were held
identical between this study's model and the `pvlib` calls for a fair
comparison isolating the transposition/geometry formulas themselves.

## Results — after both fixes (longitude + AoI; see `poa_azimuth_fix_notes.md`)

| Site | This study (`G_front`, kWh/m²/yr) | pvlib isotropic | Δ vs isotropic | pvlib Perez | Δ vs Perez |
|---|---|---|---|---|---|
| Konya | 1,939.6 | 1,939.7 | -0.00% | 2,001.2 | -3.08% |
| Almería | 2,102.6 | 2,103.4 | -0.04% | 2,170.8 | -3.14% |
| Ouagadougou | 2,259.2 | 2,259.6 | -0.02% | 2,294.5 | -1.54% |
| Freiburg | 1,458.4 | 1,459.3 | -0.06% | 1,517.9 | -3.92% |

## Interpretation

**Agreement with pvlib's isotropic model is essentially exact (≤0.06% at
every site).** This is expected, not a coincidence: after the two fixes,
this study's angle-of-incidence formula (Duffie & Beckman closed-form)
and `pvlib`'s isotropic transposition model are mathematically equivalent
formulations of the same physics. This confirms the geometric/AoI portion
of the POA model is now correct, independent of this study's own
solar-position code (since `pvlib`'s solar position comes from the
independent NREL SPA algorithm, not this study's Spencer-formula
implementation).

**The residual gap versus Perez (1.5–3.9%) is genuine model uncertainty,
not a bug.** It reflects the choice of isotropic vs. anisotropic
sky-diffuse modeling — the isotropic model treats diffuse sky radiation
as uniform across the sky dome, while Perez apportions extra weight to
circumsolar and horizon-brightening regions. This is a well-documented,
bounded source of uncertainty in the PV modeling literature and is small
relative to the ~14–16% bias the original (pre-fix) AoI formula
introduced. No further code change is recommended on this basis; it is
appropriately treated as a stated model limitation.

## Before-fix comparison (for reference)

Run with the same `pvlib` methodology but against the pre-fix
`compute_POA()`/`compute_solar_angles()` (missing longitude correction):

| Site | This study (pre-fix) | pvlib isotropic | Δ |
|---|---|---|---|
| Konya | 1,752.4 | 1,939.7 | **-9.65%** |
| Almería | 2,103.6 | 2,103.4 | +0.01% |
| Ouagadougou | 2,259.1 | 2,259.6 | -0.02% |
| Freiburg | 1,441.2 | 1,459.3 | -1.24% |

The large, site-specific discrepancy at Konya (longitude 32.49°E, the
largest of the four sites) versus near-perfect agreement at Almería and
Ouagadougou (longitude near 0°) was the diagnostic signal that led to
discovering the missing longitude correction — see
`poa_azimuth_fix_notes.md` for the full diagnosis.

## Reproducing this validation

```python
import pvlib
import numpy as np

solpos = pvlib.solarposition.get_solarposition(times_utc, lat, lon)
dni_extra = pvlib.irradiance.get_extra_radiation(times_utc)
airmass = pvlib.atmosphere.get_relative_airmass(solpos['apparent_zenith'])

common = dict(
    surface_tilt=beta, surface_azimuth=180,
    solar_zenith=solpos['apparent_zenith'].values,
    solar_azimuth=solpos['azimuth'].values,
    dni=df["DNI"].values, ghi=df["GHI"].values, dhi=df["DHI"].values,
    dni_extra=dni_extra.values, airmass=airmass.values, albedo=0.25,
)
iso = pvlib.irradiance.get_total_irradiance(model='isotropic', **common)
perez = pvlib.irradiance.get_total_irradiance(model='perez', **common)
```

`pvlib` version used: 0.15.2. Requires `DNI`, `GHI`, `DHI` from the site's
TMY CSV and the site's `(lat, lon)` from `config_00.SITES`.

## Scope and limitations of this validation

This validates the **POA irradiance geometry only** — angle of incidence,
solar position, and sky-diffuse transposition. It does not independently
validate:
- Cell temperature (Faiman model) — no external benchmark run
- PV power conversion (module efficiency, system losses) — no external
  benchmark run
- Bifacial rear-side irradiance — `pvlib` was run in front-side-only mode
  for this comparison; the bifacial gain itself (`phi = 0.70`) is a fixed
  manufacturer-specified parameter (Table 4), not separately validated
  against a bifacial-specific external tool

A full-chain validation against `pvlib.pvsystem` (including temperature
and power models) is a reasonable extension for future work but was
outside the scope of what the geometric bug discovery required.
