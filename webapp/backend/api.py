"""
Minimal FastAPI wrapper around the real, unmodified APV+AD pipeline.

Endpoints (see docs/webapp_africa_platform_design.md section 4 for the
full sketch; /api/pareto is not built yet):

  GET  /api/climate?lat=&lon=          -- live PVGIS TMY summary
  GET  /api/sites/reference            -- the 4 published sites' real results
  POST /api/evaluate                   -- single eLER_objective() call, any point
  POST /api/optimize                   -- real DE + L-BFGS-B run, as a background job
  GET  /api/optimize/{job_id}          -- poll job status/progress/result

Run locally:
  uvicorn api:app --reload --port 8000

Concurrency: the site definition for a given request is scoped to that
request's context (see model_bridge.runtime_site); concurrent requests on
the same coordinates no longer collide. JOBS itself is still a
process-global in-memory dict -- fine for a single-worker dev server, but
it must move to shared storage before running multiple workers.
"""
import copy
import sys
import threading
import uuid
from pathlib import Path
from typing import Optional

import pandas as pd
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "notebooks"))
import config_00  # noqa: E402

from model_bridge import (  # noqa: E402
    evaluate_catalog_site, evaluate_point, OUTPUTS_CSV_DIR, UnknownSiteError,
)
from pvgis_client import fetch_tmy, add_growing_season  # noqa: E402
from economics import compute_4e  # noqa: E402
from optimizer import optimize_site  # noqa: E402
from pareto import run_pareto  # noqa: E402

JOBS: dict[str, dict] = {}

app = FastAPI(title="APV+AD Atlas API", version="0.1.0")

# Local-dev only: allows webapp/frontend/index.html (opened as a file:// page,
# or served from any port) to call this API directly from the browser.
# Tighten allow_origins before this ever leaves a single developer's machine.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.exception_handler(UnknownSiteError)
def _unknown_site_handler(request, exc: UnknownSiteError):
    """
    Surface a site-resolution failure as an explicit 422, never as a bare
    KeyError escaping to an undocumented HTTP 500. Before the 2026-09
    concurrency fix this was the dominant failure mode under load -- see
    model_bridge's module docstring and audit/phase0/.
    """
    return JSONResponse(
        status_code=422,
        content={"detail": f"Site resolution failed: {exc.args[0] if exc.args else exc}"},
    )

CATALOG_SITES = list(config_00.SITES.keys())  # ["Konya", "Almeria", "Ouagadougou", "Freiburg"]

# Crop presets, matching webapp/design/Analyze.dc.html's CROP_META and
# docs/webapp_africa_platform_design.md section 3. ONLY PAR_sat is
# defensibly derived here (a standard PPFD light-saturation point, same
# conversion basis as tomato's 174 W/m^2). Residue ratio / dry matter /
# yield for non-tomato crops are NOT sourced -- do not add them here
# without a real citation; leave the crop borrowing its host site's
# values for everything else, same honesty rule as demo_new_point.py.
CROP_PRESETS = {
    "tomato":  {"PAR_sat": 174.0, "confidence_tier": "validated",         "source": "Mohammedi et al. 2023"},
    "lettuce": {"PAR_sat": 90.0,  "confidence_tier": "literature_estimated", "source": "generic C3 PPFD saturation point, not independently validated"},
    "cassava": {"PAR_sat": 150.0, "confidence_tier": "literature_estimated", "source": "generic C3 PPFD saturation point, not independently validated"},
    "maize":   {"PAR_sat": 400.0, "confidence_tier": "literature_estimated", "source": "generic C4 PPFD saturation point, not independently validated"},
}

DEFAULT_DESIGN_VARS = [16.37, 3.026, 1.501, 2.64, 31.8, 0.0502]  # Ouagadougou's optimized design, used as a generic starting point


class EvaluateRequest(BaseModel):
    lat: float
    lon: float
    crop: str = "tomato"
    scenario: str = "S0"
    borrow_from: str = "Ouagadougou"          # which catalog site's crop/substrate/econ defaults to borrow
    planting_mmdd: Optional[str] = None       # e.g. "07-01" -- defaults to borrow_from's own dates
    harvest_mmdd: Optional[str] = None
    design_vars: Optional[list[float]] = None  # [beta, d_row, H_m, V_dig, HRT, f_PV_heat]


@app.get("/api/climate")
def get_climate(lat: float = Query(...), lon: float = Query(...)):
    try:
        df = fetch_tmy(lat, lon)
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"PVGIS fetch failed: {exc}")
    return {
        "lat": lat, "lon": lon,
        "elevation_m": df.attrs["elevation_m"],
        "radiation_db": df.attrs["radiation_db"],
        "annual_ghi_kwh_m2": round(df["GHI"].sum() / 1000, 1),
        "mean_temp_c": round(df["T2m"].mean(), 1),
        "mean_wind_ms": round(df["WS"].mean(), 1),
    }


@app.get("/api/sites/reference")
def get_reference_sites(scenario: str = "S0"):
    """The 4 published sites' real Table 8/9/12 values -- trust anchors."""
    path = OUTPUTS_CSV_DIR / "results_summary.csv"
    if not path.exists():
        raise HTTPException(status_code=404, detail="results_summary.csv not found")
    df = pd.read_csv(path)
    rows = df[df["scenario"] == scenario]
    return rows.to_dict(orient="records")


@app.post("/api/evaluate")
def post_evaluate(req: EvaluateRequest):
    if req.borrow_from not in CATALOG_SITES:
        raise HTTPException(status_code=400, detail=f"borrow_from must be one of {CATALOG_SITES}")
    if req.crop not in CROP_PRESETS:
        raise HTTPException(status_code=400, detail=f"crop must be one of {list(CROP_PRESETS)}")

    host = config_00.SITES[req.borrow_from]
    planting = req.planting_mmdd or host["planting_date"]
    harvest = req.harvest_mmdd or host["harvest_date"]
    crop_preset = CROP_PRESETS[req.crop]

    try:
        df = fetch_tmy(req.lat, req.lon)
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"PVGIS fetch failed: {exc}")
    df = add_growing_season(
        df, planting, harvest,
        source=f"borrowed from {req.borrow_from}" if not (req.planting_mmdd or req.harvest_mmdd) else "user-supplied",
    )

    site_dict = copy.deepcopy(host)
    site_dict["lat"] = req.lat
    site_dict["lon"] = req.lon
    site_dict["elevation_m"] = df.attrs["elevation_m"]
    site_dict["PAR_sat"] = crop_preset["PAR_sat"]

    borrowed_fields = sorted(
        k for k in site_dict
        if k not in ("lat", "lon", "elevation_m", "PAR_sat")
    )
    if req.crop != "tomato":
        borrowed_fields = ["PAR_sat (crop-specific, literature-estimated)"] + borrowed_fields

    x = req.design_vars or DEFAULT_DESIGN_VARS
    result = evaluate_point(f"_runtime_{req.lat}_{req.lon}", site_dict, df, x, scenario=req.scenario,
                             PAR_sat=crop_preset["PAR_sat"])
    if result is None:
        raise HTTPException(status_code=422, detail="Design vector infeasible for this site (check min_beta_deg / row spacing / OLR bounds)")

    fourE = compute_4e(site_dict, result, scenario=req.scenario)

    return {
        "location": {"lat": req.lat, "lon": req.lon, "elevation_m": site_dict["elevation_m"]},
        "crop": req.crop,
        "crop_confidence_tier": crop_preset["confidence_tier"],
        "borrowed_from": req.borrow_from,
        "borrowed_fields": borrowed_fields,
        "season_source": df.attrs.get("season_source"),
        "design_vars_source": "user-supplied" if req.design_vars else f"reused from {req.borrow_from}'s optimized design (not re-optimized for this point)",
        "economics_note": f"NPV/IRR/LCOE use {req.borrow_from}'s prices, discount rate, and grid emission factor -- Bamako-specific economics would need their own sourced values, same honesty rule as the physical borrowed_fields above.",
        "result": result,
        "economics": fourE,
    }


class OptimizeRequest(BaseModel):
    lat: float
    lon: float
    crop: str = "tomato"
    scenario: str = "S0"
    borrow_from: str = "Ouagadougou"
    planting_mmdd: Optional[str] = None
    harvest_mmdd: Optional[str] = None
    max_generations: Optional[int] = None  # override DE_SETTINGS default (50) for faster demo runs


def _run_optimize_job(job_id: str, req: OptimizeRequest) -> None:
    job = JOBS[job_id]
    try:
        host = config_00.SITES[req.borrow_from]
        planting = req.planting_mmdd or host["planting_date"]
        harvest = req.harvest_mmdd or host["harvest_date"]
        crop_preset = CROP_PRESETS[req.crop]

        df = fetch_tmy(req.lat, req.lon)
        df = add_growing_season(
            df, planting, harvest,
            source=f"borrowed from {req.borrow_from}" if not (req.planting_mmdd or req.harvest_mmdd) else "user-supplied",
        )

        site_dict = copy.deepcopy(host)
        site_dict["lat"] = req.lat
        site_dict["lon"] = req.lon
        site_dict["elevation_m"] = df.attrs["elevation_m"]
        site_dict["PAR_sat"] = crop_preset["PAR_sat"]

        def progress(gen, best):
            job["progress"].append({"generation": gen, "best_eLER": round(best, 4)})

        x_opt, full, history = optimize_site(
            f"_job_{job_id}", site_dict, df, scenario=req.scenario,
            PAR_sat=crop_preset["PAR_sat"], progress_cb=progress,
            max_generations=req.max_generations,
        )
        fourE = compute_4e(site_dict, full, scenario=req.scenario)
        job["status"] = "done"
        job["result"] = {
            "location": {"lat": req.lat, "lon": req.lon, "elevation_m": site_dict["elevation_m"]},
            "crop": req.crop,
            "crop_confidence_tier": crop_preset["confidence_tier"],
            "borrowed_from": req.borrow_from,
            "design_vars_optimized": x_opt,
            "result": full,
            "economics": fourE,
        }
    except UnknownSiteError as exc:
        job["status"] = "error"
        job["error_kind"] = "configuration"
        job["error"] = str(exc)
    except Exception as exc:  # noqa: BLE001 -- surfaced to the client via job["error"]
        job["status"] = "error"
        job["error_kind"] = "internal"
        job["error"] = str(exc)


@app.post("/api/optimize")
def post_optimize(req: OptimizeRequest):
    if req.borrow_from not in CATALOG_SITES:
        raise HTTPException(status_code=400, detail=f"borrow_from must be one of {CATALOG_SITES}")
    if req.crop not in CROP_PRESETS:
        raise HTTPException(status_code=400, detail=f"crop must be one of {list(CROP_PRESETS)}")

    job_id = uuid.uuid4().hex[:12]
    JOBS[job_id] = {"status": "running", "progress": [], "result": None, "error": None}
    threading.Thread(target=_run_optimize_job, args=(job_id, req), daemon=True).start()
    return {"job_id": job_id}


@app.get("/api/optimize/{job_id}")
def get_optimize(job_id: str):
    job = JOBS.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="unknown job_id")
    return job


class ParetoRequest(BaseModel):
    lat: float
    lon: float
    crop: str = "tomato"
    borrow_from: str = "Ouagadougou"
    planting_mmdd: Optional[str] = None
    harvest_mmdd: Optional[str] = None
    pop_size: int = 50
    n_gen: int = 80


def _run_pareto_job(job_id: str, req: ParetoRequest) -> None:
    job = JOBS[job_id]
    try:
        host = config_00.SITES[req.borrow_from]
        planting = req.planting_mmdd or host["planting_date"]
        harvest = req.harvest_mmdd or host["harvest_date"]
        crop_preset = CROP_PRESETS[req.crop]

        df = fetch_tmy(req.lat, req.lon)
        df = add_growing_season(df, planting, harvest,
                                 source=f"borrowed from {req.borrow_from}" if not (req.planting_mmdd or req.harvest_mmdd) else "user-supplied")

        site_dict = copy.deepcopy(host)
        site_dict["lat"] = req.lat
        site_dict["lon"] = req.lon
        site_dict["elevation_m"] = df.attrs["elevation_m"]
        site_dict["PAR_sat"] = crop_preset["PAR_sat"]

        def progress(gen, n_gen):
            job["progress"].append({"generation": gen, "n_gen": n_gen})

        eLER, LER_water, X = run_pareto(f"_pareto_{job_id}", site_dict, df,
                                         pop_size=req.pop_size, n_gen=req.n_gen,
                                         PAR_sat=crop_preset["PAR_sat"], progress_cb=progress)
        order = eLER.argsort()
        job["status"] = "done"
        job["result"] = {
            "location": {"lat": req.lat, "lon": req.lon, "elevation_m": site_dict["elevation_m"]},
            "crop": req.crop,
            "crop_confidence_tier": crop_preset["confidence_tier"],
            "borrowed_from": req.borrow_from,
            "front": [
                {"eLER_defA": round(float(eLER[i]), 3), "LER_water_season_pct": round(float(LER_water[i]), 2)}
                for i in order
            ],
        }
    except UnknownSiteError as exc:
        job["status"] = "error"
        job["error_kind"] = "configuration"
        job["error"] = str(exc)
    except Exception as exc:  # noqa: BLE001
        job["status"] = "error"
        job["error_kind"] = "internal"
        job["error"] = str(exc)


@app.post("/api/pareto")
def post_pareto(req: ParetoRequest):
    if req.borrow_from not in CATALOG_SITES:
        raise HTTPException(status_code=400, detail=f"borrow_from must be one of {CATALOG_SITES}")
    if req.crop not in CROP_PRESETS:
        raise HTTPException(status_code=400, detail=f"crop must be one of {list(CROP_PRESETS)}")

    job_id = uuid.uuid4().hex[:12]
    JOBS[job_id] = {"status": "running", "progress": [], "result": None, "error": None}
    threading.Thread(target=_run_pareto_job, args=(job_id, req), daemon=True).start()
    return {"job_id": job_id}


@app.get("/api/pareto/{job_id}")
def get_pareto(job_id: str):
    job = JOBS.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="unknown job_id")
    return job
