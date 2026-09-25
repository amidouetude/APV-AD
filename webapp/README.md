# webapp/

The Africa APV+AD platform (design + backend), built alongside but never
modifying the research repository one level up. See
`../docs/webapp_africa_platform_design.md` for the full design brainstorm
and build order this follows.

## What's here

- **`design/`** -- a clickable design-canvas prototype (three screens: Map,
  Analyze, Compare) published as a Claude Artifact. Static mockups only;
  no connection to `backend/`. Source files are the `.dc.html` artboards +
  `canvas.json`; `apv-ad-atlas-prototype.html` is the seeded, published copy.
- **`backend/`** -- a real, working API wrapping the actual validated model:
  - `model_bridge.py` -- imports `notebooks/simulation_functions.py` and
    `notebooks/config_00.py` **unmodified**. Generalizes from "4 fixed
    sites" to "any point" by registering a throwaway entry into
    `config_00.SITES` at call time (`register_runtime_site` /
    `evaluate_point`), never editing either file on disk.
  - `pvgis_client.py` -- live PVGIS-SARAH3 TMY fetch + the DataFrame
    transform the pipeline needs (including the DatetimeIndex it requires
    -- see the docstring in `model_bridge.load_hourly_csv` for why this
    matters and what breaks silently without it).
  - `validate_ouagadougou.py` -- the correctness gate: reproduces
    `outputs/csv/results_summary.csv`'s Ouagadougou/S0 row from the live
    code path. Currently passes on every eLER/energy/water field; three
    fields (PV_sold, ESR, W_saved) differ by <0.3%, attributable to
    reusing the CSV's *rounded* design vector rather than the optimizer's
    full-precision floats.
  - `demo_new_point.py` -- runs the real pipeline on live PVGIS data for
    Bamako, Mali (never part of the study), printing which outputs are
    genuinely site-specific vs. borrowed from Ouagadougou's crop/
    substrate/economic assumptions.
  - `economics.py` -- the 4E economics (CAPEX, revenue by scenario, NPV,
    IRR, LCOE, CO2 avoidance), extracted by hand from
    `notebooks/04_metrics.ipynb` cells 2-5 (formulas unchanged, `SITES[name]`
    lookups replaced by a `site` dict parameter). Validated against
    `results_summary.csv` for Ouagadougou under both S0 and S6 -- see
    `validate_economics.py`. **Not auto-synced** with the notebook (unlike
    `simulation_functions.py`, which `02_simulation.ipynb` exports
    automatically) -- if the notebook's formulas change, this file needs
    the same edit by hand.
  - `api.py` -- FastAPI service: `GET /api/climate`,
    `GET /api/sites/reference`, `POST /api/evaluate` (now returns both the
    physical/eLER result AND the 4E economics for any point). Tested via
    `TestClient` (see below); not yet run as a standalone server, no
    `/api/optimize` (DE) or `/api/pareto` (NSGA-II) endpoint yet.

## Running the backend

```bash
pip install -r webapp/backend/requirements.txt   # fastapi, uvicorn, httpx2
cd webapp/backend
python validate_ouagadougou.py     # correctness gate
python demo_new_point.py           # live PVGIS + a genuinely new location
uvicorn api:app --reload --port 8000   # serve the API
```

## Status vs. the design doc's build order

1. ~~`site: str | dict` refactor~~ -- **not needed**: the runtime-injection
   approach in `model_bridge.py` achieves the same generalization without
   touching `simulation_functions.py`, which better satisfies "don't
   modify the existing files" than the originally proposed signature change.
2. Data ingestion service -- **done** (`pvgis_client.py`), live-tested
   against Bamako.
3. `/api/evaluate` fast tier -- **done**, tested for both a `validated`
   crop (tomato) and a `literature_estimated` one (maize), correctly
   surfacing the confidence tier and the list of borrowed (non-sourced)
   fields in the response.
4. `/api/optimize` (DE job + live convergence) -- **not started**.
5. Map layer wired to the backend -- **not started** (the `design/`
   prototype's map is still static sample data).
6. Multi-crop / multi-substrate presets -- **partially started**: `api.py`
   has a small `CROP_PRESETS` dict (tomato validated, lettuce/cassava/maize
   literature-estimated PAR_sat only -- residue ratio/DM/yield for those
   three are deliberately NOT filled in, per the honesty rule in section 6
   of the design doc). No substrate presets yet.
7. `/compare` page wired to real data -- **not started**.

## Resolved gap

~~The 4E economics exist only inside `04_metrics.ipynb`~~ -- extracted into
`economics.py` and validated against Ouagadougou's published S0 and S6 rows.
`/api/evaluate` now returns eLER + energy + water + full 4E economics for
any point (see `demo_new_point.py`-style output above for Bamako, Mali:
eLER 1.736, NPV $788k, IRR 27.8% (S0) / 29.6% (S6) -- economics honestly
flagged as using Ouagadougou's borrowed prices/discount rate/grid EF, not
Mali-specific ones).

## Resolved gap: real optimizer

`optimizer.py` -- the real `scipy.optimize.differential_evolution` +
L-BFGS-B pipeline from `notebooks/03_optimizer.ipynb` cells 2-3 (identical
hyperparameters from `config_00.DE_SETTINGS`), generalized to any site via
the same runtime-injection mechanism as `model_bridge.py`. Validated
(`validate_optimizer.py`) against Ouagadougou: converges to **eLER 1.7230
vs. published 1.7230**, every design variable matching to 3-4 decimals,
including the paper's own boundary-convergence finding (row spacing and
module height landing on their physical minimums). ~200s runtime, matching
the paper's reported 102-172s.

Exposed via `api.py` as an async job: `POST /api/optimize` returns a
`job_id` immediately; `GET /api/optimize/{job_id}` polls status, a live
per-generation progress list, and the final optimized result + economics
once done. Tested end-to-end with a short `max_generations` override for
a fast sanity check of the wiring (job creation -> background thread ->
polling -> completion) -- the full 50-generation convergence quality is
what `validate_optimizer.py` already confirmed separately.

**Known limitation**: `JOBS` and the site-injection it relies on are
process-global in-memory state -- fine for a single-user dev server, not
safe for concurrent multi-worker deployment as-is (same caveat already
noted on `model_bridge.register_runtime_site`).

## Resolved gap: a real (unsandboxed) frontend

**Important constraint discovered here**: the `design/` canvas prototype
published earlier is a Claude Design Artifact, which runs in a sandboxed
iframe with **no network egress** except its own origin and Google Fonts
-- a hard platform rule, not a bug. It can never call `backend/api.py`,
however it's wired. "Connecting the prototype to the real API" therefore
has to mean a separate, genuinely unsandboxed page.

- **`frontend/index.html`** -- a real, working page (no build step, no
  framework) that calls the live backend directly: location presets
  (the 4 published sites + two brand-new cities, Bamako and Accra), a
  crop selector, "Get live climate" (`/api/climate`), "Evaluate"
  (`/api/evaluate`), and "Run full optimization" (`/api/optimize` +
  polling, with a live generation counter while the real ~2-4 min DE run
  is in progress). Visually modeled on the `design/` prototype's
  Academic Paper Light theme, but every number on this page is real.
- `api.py` now has `CORSMiddleware` (`allow_origins=["*"]`) so a
  `file://`-opened page can call `localhost:8000` -- verified via a raw
  CORS preflight check, not just assumed.
- **Not verified in an actual browser** (none available in this
  environment) -- the full data/network contract (CORS, request/response
  shapes) was verified via `curl`/`TestClient` against the exact calls
  `index.html` makes, but the page's on-screen rendering hasn't been
  visually confirmed. Open `frontend/index.html` directly to check.

## Running the full live stack

```bash
pip install -r webapp/backend/requirements.txt
cd webapp/backend
uvicorn api:app --host 127.0.0.1 --port 8000   # leave running
```
Then open `webapp/frontend/index.html` directly in a browser (no server
needed for the frontend itself -- plain `file://` works because of CORS).

## Resolved gap: a real map

`frontend/index.html` now embeds a genuine Leaflet + OpenStreetMap map of
Africa (not the schematic silhouette from `design/`, since this page is
unsandboxed and can load real tiles) above the existing evaluate/optimize
form:

- Diamond markers for the 4 published reference sites, popups showing
  their real eLER/IRR/NPV from `/api/sites/reference`.
- Click anywhere in Africa -> drops a pin, calls `/api/climate` live, and
  shows real PVGIS annual GHI/temp/wind in the popup -- this substitutes
  for a precomputed choropleth layer (which would need its own raster
  data pipeline, not built) with an on-demand live lookup instead.
- "Use this location" on any popup (reference or dropped pin) fills the
  lat/lon fields, so map -> evaluate/optimize is now one continuous flow
  in a single real page, closing the original "map + analyze" loop this
  whole project started from.

**Visually verified**, not just structurally: rendered headlessly via
Chrome (`chrome.exe --headless=new --screenshot`, no Playwright/Selenium
needed) against the live backend. First render caught a real bug --
the hardcoded initial `setView` clipped the Freiburg marker (47.99N)
off the top of the 400px map at low zoom -- fixed by `fitBounds()` on
the actual reference-site coordinates instead. A second render (with an
auto-click of "Evaluate" injected into a throwaway copy, since Chrome's
one-shot screenshot mode can't script interactions) confirmed the full
result card renders correctly: eLER 1.736, NPV $788.1k, IRR 27.77%,
stacked bar, tiles, and the borrowed-fields honesty note all matching
the API contract exactly.

## Resolved gap: the NSGA-II Pareto front

`pareto.py` -- the real bi-objective NSGA-II search (maximise eLER and
`LER_water` jointly) from `notebooks/07_pareto_analysis.py`, generalized
to any site via the same runtime-injection mechanism. Validated
(`validate_pareto.py`) against Konya: front extremes come out at
**eLER 1.202/1.641 vs. published 1.20/1.64**, **LER_water 36.1%/18.9% vs.
published 36.1%/18.9%** (paper Section 5.3) -- both ends of the trade-off
curve reproduce exactly, including that the eLER-maximising end matches
the scalar DE optimum, as the paper itself notes. ~250s runtime for the
full 50-pop/80-gen search (4,000 evaluations).

Exposed as `POST /api/pareto` (job) / `GET /api/pareto/{job_id}` (poll),
same async pattern as `/api/optimize`. Tested end-to-end with a tiny
pop/gen override for a fast wiring check. Not yet surfaced in
`frontend/index.html` -- no Pareto-front chart in the UI yet.

## Next gaps

- No country-level economic defaults table -- every non-catalog point
  still must borrow a full site's economics via `borrow_from`. Deliberately
  not built yet: filling it in would mean asserting country-specific
  electricity prices/grid emission factors I can't verify or cite from
  this environment, exactly the kind of unsourced-data risk flagged
  throughout this project's own design docs.
- No precomputed wind/solar choropleth layer (the live-click-to-fetch
  pattern already covers the same need without one).
- `frontend/index.html` has no UI for the new `/api/pareto` endpoint yet
  (map, evaluate, and optimize are wired; Pareto is backend-only so far).
