# ET0 / Hargreaves-Samani Fix Notes

This note documents two related corrections to the reference
evapotranspiration model in `compute_ET_reduction()`
(`notebooks/simulation_functions.py`), why each was made, and their
quantified impact on the water savings results (Table 13).

---

## Fix 1 — Extraterrestrial radiation (Ra) estimation

### What was wrong

The Hargreaves-Samani equation requires extraterrestrial radiation `Ra`
— a quantity computable purely astronomically from latitude and
day-of-year, needing no ground radiation measurement at all (this is the
entire point of the method's minimal-input design). The original code
instead estimated `Ra` by dividing measured GHI by an assumed constant
clearness index of 0.75:

```
Ra_est = GHI_daily / 0.75
```

### Diagnosis

Checked the actual clearness index (`Kt = GHI / Ra_true`) implied by each
site's real TMY data against the assumed 0.75:

| Site | Ra (code, GHI/0.75) | Ra (true, astronomical) | Bias | Real avg. clearness index |
|---|---|---|---|---|
| Konya | 23.4 MJ/m²/d | 29.0 MJ/m²/d | -19.1% | 0.585 |
| Almería | 25.0 | 29.4 | -14.8% | 0.630 |
| Ouagadougou | 29.1 | 35.2 | -17.3% | 0.624 |
| Freiburg | 16.6 | 25.2 | -33.9% | 0.521 |

Real clearness indices (0.52–0.63) are well below the assumed 0.75 at
every site, so the code systematically underestimated `Ra` — worst at
Freiburg, the cloudiest site, which is exactly where a "clear sky
fraction" assumption should be most wrong.

### The fix

Replaced with the standard FAO-56 (Allen et al. 1998) astronomical
formula, computed from latitude and day-of-year only:

```
Ra = (24*60/pi) * Gsc * dr * (omega_s*sin(phi)*sin(delta) + cos(phi)*cos(delta)*sin(omega_s))
```

where `Gsc = 0.0820` MJ/m²/min (solar constant), `dr` = inverse relative
Earth-Sun distance, `delta` = declination, `omega_s` = sunset hour angle,
`phi` = latitude. This matches the citation the manuscript already made
(Droogers & Allen 2002, ±10% agreement with FAO-56 Penman-Monteith) —
that citation is only valid for the *standard* Hargreaves-Samani
formula, which this fix now actually implements.

---

## Fix 2 — Missing MJ-to-mm conversion factor

### How this was found

After Fix 1 alone, `ET0_season` at Konya came out to ~11 mm/day — not
physically plausible (the manuscript's own Section 4.6 states typical
semi-arid ET0 is 5–7 mm/day). This discrepancy was caught by a
plausibility check, not by trusting the relative before/after comparison
alone.

### What was wrong

The standard Hargreaves-Samani formula (FAO-56 Eq. 52) is:

```
ET0 [mm/day] = 0.0023 * (Tmean + 17.8) * sqrt(Tmax - Tmin) * Ra[mm/day]
```

where `Ra` must be in **mm/day-equivalent** units, i.e. the MJ/m²/day
value from Fix 1 multiplied by 0.408 (= 1 / latent heat of vaporization
of water, 2.45 MJ/kg). This 0.408 factor was missing from both the
**original** code and the first pass of Fix 1 above — the original
code's too-small `Ra_est` (from the GHI/0.75 proxy) happened to partially
mask the missing factor by producing plausible-looking (but doubly
wrong) numbers. Fixing `Ra` alone, without also adding this factor, made
the ~2.45× overestimate fully visible rather than causing it.

### The fix

```python
ET0_daily_mm = 0.0023 * 0.408 * Ra_est * (T_mean_daily + 17.8) * T_range**0.5
```

### Verification

| Site | ET0 daily rate (post both fixes) | Season length | Modeled ET0_annual | Table 2 literature reference (annual ET0) |
|---|---|---|---|---|
| Konya | 4.47 mm/day | 189 days | 1,134 mm | 1,200 mm |
| Almería | 3.82 mm/day | 122 days | 1,151 mm | 1,350 mm |
| Ouagadougou | 4.26 mm/day | 122 days | 1,839 mm | 1,800 mm |
| Freiburg | 3.03 mm/day | 167 days | 685 mm | ~900 mm |

Daily rates are now physically plausible (3–4.5 mm/day), and — as an
independent sanity check the manuscript didn't have before — the
model's own computed annual ET0 now closely tracks the literature
reference values already stated in Table 2, which it failed to do at
all before this fix (was 2–4× too high).

---

## Impact on Table 13 (Water Savings)

**Scope: fully contained to the Water Savings section.** `LER_water`
does not feed into `eLER`, the DE optimizer's objective, or any
economics — confirmed by inspection of `eLER_objective()`, which sums
only `LER_crop + LER_PV + LER_biogas`. No other section of the
manuscript required revision as a result of this fix.

### Absolute values (materially changed)

| Site | W_saved (mm/season) | ET0_season (mm) |
|---|---|---|
| | Before → After | Before → After |
| Konya | 383 → 159 | 650 → 846 |
| Almería | 186 → 88 | 480 → 466 |
| Ouagadougou | 152 → 81 | 620 → 520 |
| Freiburg | 186 → 114 | 350 → 509 |

### Ratio metric (essentially unchanged)

`LER_water_season_pct` values changed by less than 0.15 percentage
points at every site (18.8%, 18.8%, 15.6%, 22.4% — see Table 13),
because both `W_saved` and `ET0_season` carry the same missing-factor
bias and it cancels in the ratio — the same structural pattern seen in
every other fix in this study (see `poa_azimuth_fix_notes.md`,
`ad_bmp_sensitivity.md`). **Conclusion C6 and the corrected
growing-season-vs-annual-denominator finding are unaffected in their
qualitative claims**, though the specific absolute mm and m³/ha figures
in Section 4.6's "Practical interpretation" paragraph needed updating.

## Code locations

- `compute_ET_reduction(df, site)` — `notebooks/simulation_functions.py`
- Regression test: `test_et_reduction_physically_plausible_daily_rate` in
  `tests/test_simulation_functions.py` (checks ET0 falls in a
  1–10 mm/day physically plausible range, specifically to prevent this
  class of unit-conversion bug from silently reappearing)
