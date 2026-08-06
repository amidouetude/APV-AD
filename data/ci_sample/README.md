# CI Sample Data

## `hourly_Konya_sample.csv`

A 720-hour (April 2005) slice of the real Konya PVGIS-SARAH3 TMY data used
throughout this study — not synthetic data. Extracted from the full
8760-hour `outputs/csv/hourly_Konya.csv` produced by `01_load_data.ipynb`
by filtering to `month == 4`.

**Why this slice:** it spans Konya's April 25 planting-date boundary,
giving 144 in-season hours and 576 off-season hours in a single small
(~60 KB) file — enough to exercise both code branches of the
season-dependent logic (crop yield, ET reduction, AD residue input) with
genuine weather data, without committing the full multi-hundred-KB
per-site, 8760-hour dataset.

**Used by:** `tests/test_notebook_smoke.py`, for a real-data end-to-end
pipeline smoke test in CI (see that file's module docstring for the
rationale).

## Regenerating this file

```python
import pandas as pd
df = pd.read_csv("outputs/csv/hourly_Konya.csv", index_col=0, parse_dates=True)
april = df[df["month"] == 4]
april.to_csv("data/ci_sample/hourly_Konya_sample.csv")
```

## Data source and license

Original data: PVGIS-SARAH3 Typical Meteorological Year database,
European Commission Joint Research Centre (JRC), accessed via
[PVGIS](https://re.jrc.ec.europa.eu/pvg_tools/en/). PVGIS data is made
available by the JRC under a permissive reuse policy requiring
attribution (see the [PVGIS terms of
use](https://joint-research-centre.ec.europa.eu/photovoltaic-geographical-information-system-pvgis/getting-started-pvgis/pvgis-releases/pvgis-53_en)
for current terms). This sample is a small derived excerpt retained for
software-testing purposes; for the full dataset and citation details, see
the manuscript's Section 3.1 and `huld2012` / `pvgis` references.
