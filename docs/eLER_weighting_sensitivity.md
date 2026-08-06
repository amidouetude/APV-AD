# eLER Component Weighting: Justification and Sensitivity Analysis

This note documents why eLER uses equal weighting of its three
components by default, and the sensitivity analysis performed to test
that choice against alternative schemes, in response to peer review item
#1 (eLER sums heterogeneous ratios with no stated weighting
justification).

## The critique

```
eLER = LER_crop + LER_PV + LER_biogas
```

is a literal, unweighted sum of three ratios that are dimensionless but
represent different physical quantities: a PAR-integral ratio, a PV
energy ratio, and a biogas energy ratio. The reviewer's concern: treating
a 1% gain in each as equally valuable to the objective has no stated
justification.

## Why equal weighting is the methodologically consistent choice

Each of the three components is already a **ratio to its own reference
system** — `LER_crop` is referenced against an open-field monoculture,
`LER_PV` against a dense-pack reference array, `LER_biogas` against a
standalone reference AD unit. This normalization is precisely what the
original two-component LER (Mead & Willey, 1980) was designed to do: by
expressing each output as a fraction of what a dedicated monoculture
would achieve on the same land, outputs that are physically
incommensurable in raw units (kg of crop, kWh of electricity, m³ of
biogas) become commensurable as land-equivalent fractions. Re-weighting
by the *raw physical magnitude* of each output re-introduces exactly the
scale problem this normalization exists to solve.

This is not just an assertion — it is demonstrated empirically below.

## Method

Two alternative weighting schemes were tested, both grounded in the
manuscript's own data (no invented parameters):

**Economic (revenue-share) weighting** — weights proportional to each
component's share of Scenario S0 gross revenue at each site:

```python
w_crop = 3 * rev_tomato / total_revenue
w_PV   = 3 * rev_pv / total_revenue
w_bio  = 3 * rev_biogas / total_revenue
```

**Energy-content weighting** — weights proportional to each component's
actual energy content (not market price), converting crop yield to
caloric energy equivalent (fresh tomato ≈ 3.10 MJ/kg = 0.861 kWh/kg,
USDA FoodData Central) and comparing directly against PV and biogas
output in MWh:

```python
w_crop = 3 * E_crop_MWh / total_energy
w_PV   = 3 * E_PV_MWh / total_energy
w_bio  = 3 * E_bio_MWh / total_energy
```

Both use `eLER_weighted()`'s normalization (weights rescaled to sum to
3, so equal weights `(1,1,1)` reproduce the unweighted default exactly).

## Results

### Component weights implied by each scheme

| Site | Scheme | w_crop | w_PV | w_bio |
|---|---|---|---|---|
| Konya | Economic | 0.266 | 2.722 | 0.012 |
| | Energy | 0.050 | 2.941 | 0.010 |
| Almería | Economic | 0.206 | 2.786 | 0.008 |
| | Energy | 0.042 | 2.949 | 0.009 |
| Ouagadougou | Economic | 0.105 | 2.886 | 0.008 |
| | Energy | 0.044 | 2.943 | 0.013 |
| Freiburg | Economic | 0.132 | 2.860 | 0.009 |
| | Energy | 0.030 | 2.961 | 0.009 |

**Both magnitude-based schemes collapse to an almost pure PV-only
metric.** PV output outweighs crop and biogas output by 50–100× in raw
physical or monetary magnitude at every site (e.g., Freiburg: PV = 938
MWh/ha vs. crop energy content = 9.6 MWh/ha). This is a physical fact
about solar panel energy density vs. photosynthesis, not a statement
about land-use value — but it means any weighting scheme built on raw
magnitude inevitably assigns PV a weight near 3 (i.e., "ignore crop and
biogas") regardless of the specific magnitude source (dollars or joules).

### Effect on site ranking (at fixed, already-optimized design)

| Weighting scheme | Ranking (highest eLER first) |
|---|---|
| Equal (default) | Freiburg > Ouagadougou > Almería > Konya |
| Energy-content | Ouagadougou > Almería > Konya > **Freiburg (last)** |
| Economic | Ouagadougou > Almería > Freiburg > Konya |

Under both magnitude-based schemes, Ouagadougou displaces Freiburg as
the top-ranked site, and under energy weighting Freiburg falls to last.

### Effect on optimal design (full DE re-optimization under economic weights)

| Site | d_row* (economic wt.) | H_m* (economic wt.) | β* (economic wt.) | d_row*/H_m*/β* (equal wt., for comparison) |
|---|---|---|---|---|
| Konya | 3.03 | 1.50 | 20.9° | 3.03 / 1.50 / 20.6° |
| Almería | 3.03 | 1.50 | 20.8° | 3.03 / 1.51 / 20.5° |
| Ouagadougou | 3.02 | 1.50 | 17.1° | 3.03 / 1.50 / 16.4° |
| Freiburg | 3.01 | 1.53 | 23.3° | 3.03 / 1.50 / 23.1° |

**Boundary convergence (Finding C3 / Section 4.1) is robust to
weighting** — `d_row` and `H_m` converge to essentially the same physical
boundary regardless of weighting scheme, since PV revenue dominates
under any economically-motivated weighting, and "maximize PV density" is
close to optimal for that objective regardless of exactly how much
weight the smaller crop/biogas terms get.

## Conclusion

1. **Design recommendations (boundary convergence, tilt angles) are
   robust** to the choice between equal and magnitude-based weighting.
2. **Cross-site eLER ranking is NOT robust** to magnitude-based
   weighting — but this instability is a property of using an
   inappropriate weighting scheme (one that reintroduces the very scale
   problem LER normalization solves), not evidence that the underlying
   metric or equal-weighting choice is arbitrary.
3. Equal weighting is retained as the primary, default metric
   (`eLER_weighted(w_crop=1, w_PV=1, w_bio=1)`, matching
   `eLER_objective()`'s default), consistent with the original LER
   framework's own normalization logic.

## Recommendation for manuscript text

Add to Section 5.3 (eLER metric interpretation) or as a new subsection:
a brief statement of the equal-weighting rationale (already added to
Methodology Section 3.8), plus a compact version of the ranking-
sensitivity table above as a robustness check, framed as: *"eLER's site
ranking is sensitive only to weighting schemes that reintroduce the
magnitude-scale problem LER normalization is designed to solve; the
optimizer's structural design recommendations (Finding C3) are unaffected
by weighting choice."*

## Code locations

- `eLER_weighted(LER_crop, LER_PV_A, LER_biogas, w_crop, w_PV, w_bio)` —
  `notebooks/simulation_functions.py`
- `eLER_objective(..., w_crop=1.0, w_PV=1.0, w_bio=1.0)` — optional
  weighting threaded through to the objective
- Regression tests: `test_eler_weighted_equal_weights_reproduces_sum`,
  `test_eler_weighted_normalizes_arbitrary_scale`,
  `test_eler_weighted_rejects_nonpositive_weight_sum` in
  `tests/test_simulation_functions.py`
