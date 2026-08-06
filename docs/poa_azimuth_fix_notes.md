# POA Angle-of-Incidence and Solar-Time Fix Notes

This note documents two related corrections to the plane-of-array (POA)
irradiance model in `compute_solar_angles()` and `compute_POA()`
(`notebooks/simulation_functions.py`), why each was made, and their
quantified impact on downstream results. Both were found and fixed during
a peer-review response cycle; the second (longitude/equation-of-time) was
discovered as a consequence of validating the first fix against `pvlib`.

---

## Fix 1 — Angle-of-incidence formula (`compute_POA`)

### What was wrong

The original angle-of-incidence formula for the tilted, south-facing PV
surface was:

```
cos(theta_i) = sin(alpha_s)*cos(beta) + cos(alpha_s)*sin(beta)
```

Algebraically this is `sin(alpha_s + beta)` (since
`sin(a)cos(b) + cos(a)sin(b) = sin(a+b)`). This expression is only exact
at solar noon, when the sun is due south (`omega = 0`). At every other
hour it silently assumes the sun is still due south, which systematically
**overestimates** `cos(theta_i)` — and therefore front-surface beam
irradiance — throughout the morning and afternoon.

### The fix

Replaced with the standard closed-form expression for a fixed, due-south
collector (Duffie & Beckman, *Solar Engineering of Thermal Processes*):

```
cos(theta_i) = sin(delta)*sin(phi - beta) + cos(delta)*cos(omega)*cos(phi - beta)
```

where `delta` = declination, `phi` = latitude, `omega` = hour angle,
`beta` = tilt. This is exact for a due-south fixed-tilt surface at every
hour of the day and requires no extra inputs beyond what
`compute_solar_angles()` already produces.

### Quantified impact (geometric, before the longitude fix below)

Computed by comparing the old and new `cos(theta_i)` formulas across every
daylight hour of a full year, at each site's actual latitude and optimal
tilt angle (Table 8):

| Site | Unweighted mean bias | Beam-weighted bias |
|---|---|---|
| Konya | +27.6% | +13.9% |
| Almería | +28.6% | +14.5% |
| Ouagadougou | +26.0% | +15.8% |
| Freiburg | +32.3% | +14.9% |

The old formula overestimated beam-weighted `cos(theta_i)` by roughly
14–16% at every site.

---

## Fix 2 — Longitude / equation-of-time correction (`compute_solar_angles`)

### How this was found

While validating Fix 1 against `pvlib`'s NREL-SPA-based solar position
(see `pvlib_validation.md`), front-surface POA matched `pvlib`'s isotropic
model to within 1.5% at three of four sites — but Konya was off by
**-9.65%**. That inconsistency, rather than the raw magnitude, was the
signal something else was wrong.

### What was wrong

`compute_solar_angles()` computed hour angle as:

```
hour_angle = (UTC_hour - 12) * 15
```

This implicitly assumes solar noon occurs at 12:00 UTC everywhere — true
only at longitude 0°E. The function never used `lon_deg` at all.

### Diagnosis

The expected hour-angle error from this omission is approximately equal
to the site's longitude in degrees (since each 15° of longitude shifts
true solar noon by 1 hour relative to UTC = 15° of hour angle):

| Site | Longitude | Expected error (deg) | Measured POA bias (pre-fix, isotropic) |
|---|---|---|---|
| Konya | 32.49°E | 32.49 | -9.65% |
| Almería | -2.46°W | -2.46 | +0.01% |
| Ouagadougou | -1.53°W | -1.53 | -0.02% |
| Freiburg | 7.85°E | 7.85 | -1.24% |

The correlation between longitude magnitude and POA bias is unmistakable
and confirms the diagnosis.

### The fix

Added longitude and equation-of-time correction to convert UTC clock time
to true local solar time before computing hour angle:

```
solar_time = hour_utc + lon_deg/15.0 + eot_min/60.0
hour_angle = (solar_time - 12) * 15
```

using the standard Spencer (1971) equation-of-time approximation.
`compute_solar_angles()` now requires `lon_deg` as a mandatory argument
(previously 2-arg, now 3-arg — a deliberate breaking change so the bug
cannot silently reproduce itself).

### Verification

After the fix, front-surface POA matches `pvlib`'s isotropic model to
within **0.06% at every site** (see `pvlib_validation.md` for the full
table) — essentially exact agreement, as expected since both formulas
are now mathematically equivalent closed-form expressions.

---

## Net impact on manuscript results (both fixes combined)

Re-running the full Differential Evolution optimizer (same settings as
Section 3.9: `NP=90`, `Gmax=50`, seed=42) with both fixes applied:

| Metric | Effect |
|---|---|
| eLER (Definition A) | Changes by ≤0.6% at every site — within DE-seed noise |
| Boundary convergence (`d_row`, `H_m`) | Unchanged — still converges to the physical constraint boundary at all sites |
| Tilt angles (`beta*`) | Shift by <1° at every site |
| Absolute PV yield | Konya +10.4% (net, after both fixes combine); other sites shift -2% to +2% |
| LCOE / NPV / IRR | Material shifts — see Results Table 12; Konya's economics improve notably since the longitude fix partially offset the AoI fix's yield reduction there |

**Why eLER is robust but absolute yield is not:** `LER_PV` is a *ratio* of
APV array output to a reference dense-pack array output, and both use the
identical (now-corrected) `compute_POA()`. The systematic bias present in
the old formula affected numerator and denominator similarly, so it
largely cancelled in the ratio. Absolute yield, CAPEX-normalized metrics
(LCOE), and cash-flow-based metrics (NPV, IRR) are not ratios of this
kind and were not protected from the bias — hence the larger shifts
there.

## Code locations

- `compute_solar_angles(df, lat_deg, lon_deg)` — `notebooks/simulation_functions.py`
- `compute_POA(df, beta_deg, lat_deg)` — `notebooks/simulation_functions.py`
- Regression tests: `test_poa_requires_lat_deg`, `test_solar_angles_requires_lon_deg`,
  `test_solar_noon_shifts_with_longitude` in `tests/test_simulation_functions.py`
