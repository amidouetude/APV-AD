# AD / BMP Parameter Sensitivity Analysis

This note documents the Monte Carlo sensitivity analysis of the anaerobic
digestion model's key parameters, performed in response to peer review
item #4 (no sensitivity analysis on `BMP_ref`, capture efficiency, or
volatile-solids fractions).

## Motivation

The AD model (`compute_AD()`, Section 3.7) uses fixed point values for
several parameters that carry real literature-reported uncertainty:
`BMP_ref` (reference biochemical methane potential), `eta_cap` (methane
capture efficiency), `VS_fraction` (crop-residue volatile-solids
fraction), and `manure_VS_fraction`. The reviewer's concern was that
reporting eLER, biogas yield, and downstream economics as point estimates
without characterizing sensitivity to these parameters overstates
precision.

## Method

`compute_AD()` and `eLER_objective()` were extended with an optional
`AD_overrides` dict (backward-compatible; `None` reproduces prior
behavior exactly — see `compute_AD`'s docstring in
`notebooks/simulation_functions.py`). This enables sampling AD parameters
without mutating global config state.

500 Monte Carlo draws per site, at each site's optimized design point
(`x*` from Table 8), sampling:

| Parameter | Distribution | Literature basis |
|---|---|---|
| `BMP_ref` | Uniform(240, 350) NmL CH₄/gVS | Manuscript §3.2: 200–280 NmL/gVS (cattle manure, Lallement et al.), 250–350 NmL/gVS (tomato residues, Lallement et al.); combined feedstock reference value 300 |
| `eta_cap` | Uniform(0.80, 1.00) | S0 baseline = 1.00 (idealized); S2 scenario tests 0.85 as "sub-optimal capture" |
| `VS_fraction` (residues) | Uniform(0.75, 0.90) | Manuscript fixes at 0.85 (Lallement et al.) |
| `manure_VS_fraction` | Uniform(0.70, 0.85) | Manuscript fixes at 0.80 (Pilarski & Pilarska) |

All four parameters were sampled simultaneously and independently
(no assumed correlation between them) per draw.

## Results

| Site | Biogas energy 90% CI (relative width) | eLER 90% CI (relative width) | eLER 90% CI (absolute width) | OLR constraint violations |
|---|---|---|---|---|
| Konya | ±44.4% | ±2.36% | 0.039 | 0 / 500 |
| Almería | ±44.9% | ±2.35% | 0.039 | 0 / 500 |
| Ouagadougou | ±43.4% | ±2.09% | 0.036 | 0 / 500 |
| Freiburg | ±43.7% | ±2.23% | 0.039 | 0 / 500 |

(90% CI computed as the p5–p95 percentile range across the 500 draws,
relative width expressed as a fraction of the median.)

## Interpretation

**Absolute biogas energy is genuinely uncertain** — a ±43–45% relative
90% CI confirms the reviewer's underlying concern was legitimate. Point
estimates of `E_biogas` (Table 10) should not be read as precise to more
than one significant figure given literature parameter uncertainty in
`BMP_ref`, `eta_cap`, and the VS fractions.

**eLER itself is robust to this uncertainty** (±2.1–2.4% relative, an
absolute range of 0.036–0.039 eLER units). This follows the same
structural pattern documented for the POA fix (see
`poa_azimuth_fix_notes.md`): `LER_biogas` is a ratio of the APV-system's
biogas output to a reference standalone AD unit's output, and both use
the same sampled AD parameters in each Monte Carlo draw. The uncertainty
in the numerator and denominator is therefore highly correlated within
a draw and cancels substantially in the ratio.

**No feasibility violations across 2,000 total draws.** Every optimized
digester design (`V_dig`, `HRT` from Table 8) remains within the OLR
operating range `[1.5, 5.0]` kg VS/m³/day even under the full literature-
cited parameter uncertainty. The optimized AD designs are not brittle
to this class of parameter uncertainty.

## Recommendation for manuscript text

Given eLER's robustness, the headline results (Tables 8–9, Section 4.2,
Conclusions C1–C3) do not need revision on the basis of this sensitivity
analysis. However, Table 10's absolute biogas energy figures, and any
economic figures downstream of biogas revenue (Scenarios S1–S6, Table
12), should be captioned or footnoted with the ±43–45% uncertainty band
to avoid over-precision. A suggested addition to Section 5 (Discussion)
or as a new Limitations item:

> AD parameter uncertainty (`BMP_ref`, capture efficiency, VS fractions)
> propagates to a ±43–45% 90% CI on absolute biogas energy, but eLER is
> robust to this uncertainty (±2–2.4%) because `LER_biogas` is a ratio
> referenced against a standalone AD unit using the same parameters,
> which largely cancels the shared uncertainty. All optimized digester
> designs remain within the OLR feasibility range across the full
> Monte Carlo sample (0/2,000 violations).

## Reproducing this analysis

```python
import numpy as np

rng = np.random.default_rng(42)
N = 500

for i in range(N):
    overrides = dict(
        BMP_ref=rng.uniform(240, 350),
        eta_cap=rng.uniform(0.80, 1.00),
        VS_fraction=rng.uniform(0.75, 0.90),
        manure_VS_fraction=rng.uniform(0.70, 0.85),
    )
    result = eLER_objective(x_star, df, site, return_full=True,
                            AD_overrides=overrides)
    # collect result["biogas_total_MWh_ha"], result["eLER_defA"],
    # result["OLR_kgVS_m3d"], etc.
```

## Code locations

- `compute_AD(df, site, V_dig_m3, HRT_days, LER_crop, scenario, AD_overrides)` —
  `notebooks/simulation_functions.py`
- `eLER_objective(..., AD_overrides=None)` — threads overrides to both the
  design-case and reference-case `compute_AD()` calls
- Regression tests: `test_compute_ad_overrides_default_matches_no_override`,
  `test_compute_ad_overrides_changes_output`,
  `test_compute_ad_higher_capture_efficiency_increases_biogas` in
  `tests/test_simulation_functions.py`
