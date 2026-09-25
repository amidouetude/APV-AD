"""
PHASE 1 -- verification des criteres d'acceptation.

1. Nord identique a 1e-9 contre les references gelees en phase 0
2. Sud coherent avec pvlib
3. lat = 0 correct
4. Zero contamination entre requetes
5. Aucune dependance MATLAB
6. Erreurs de configuration explicites
7. simulation_functions.py en phase avec le notebook (drift CI)
"""
import json, sys, threading, copy, importlib.util
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parent / "repo"
sys.path.insert(0, str(REPO / "notebooks"))
sys.path.insert(0, str(REPO / "webapp" / "backend"))
import config_00
from model_bridge import evaluate_catalog_site, evaluate_point, UnknownSiteError
from economics import compute_4e

FROZEN = Path("/mnt/user-data/outputs/phase0/c_reference_frozen.json")
ok_all = True
def check(label, ok, detail=""):
    global ok_all
    ok_all &= ok
    print(f"  [{'OK ' if ok else 'FAIL'}] {label}{(' -- ' + detail) if detail else ''}")

print("="*78); print("CRITERE 1 -- Nord identique a 1e-9 (references gelees phase 0)"); print("="*78)
frozen = json.loads(FROZEN.read_text())
worst, worst_field = 0.0, None
for site, ref in frozen.items():
    res = evaluate_catalog_site(site, ref["design_vars"], scenario="S0")
    econ = compute_4e(config_00.SITES[site], res, scenario="S0")
    for field, exp in ref["result"].items():
        if not isinstance(exp, (int, float)) or isinstance(exp, bool): continue
        got = res.get(field)
        if not isinstance(got, (int, float)): continue
        d = abs(got - exp) / max(abs(exp), 1e-12)
        if d > worst: worst, worst_field = d, f"{site}.{field}"
    for field, exp in ref["economics"].items():
        if not isinstance(exp, (int, float)) or isinstance(exp, bool): continue
        got = econ.get(field)
        if not isinstance(got, (int, float)): continue
        d = abs(got - exp) / max(abs(exp), 1e-12)
        if d > worst: worst, worst_field = d, f"{site}.econ.{field}"
check(f"ecart relatif max = {worst:.3e}", worst < 1e-9, f"pire champ: {worst_field}")

print(); print("="*78); print("CRITERE 2 -- Sud coherent avec pvlib"); print("="*78)
import pvlib
from simulation_functions import compute_POA
def cos_aoi(lat, beta, dec, ha):
    df = pd.DataFrame({"declination_deg":[float(dec)], "hour_angle_deg":[float(ha)],
                       "GHI":[0.0], "DNI":[1000.0], "DHI":[0.0], "F_shad":[0.0]})
    return float(compute_POA(df, beta, lat)["G_front"].iloc[0]) / 1000.0
rng = np.random.default_rng(11); worst_pv = 0.0; n = 0
for _ in range(2000):
    lat = float(rng.uniform(-38, -1)); beta = float(rng.uniform(0, 40))
    dec = float(rng.uniform(-23.45, 23.45)); ha = float(rng.uniform(-75, 75))
    phi, d_, w = map(np.radians, (lat, dec, ha))
    cz = np.sin(phi)*np.sin(d_) + np.cos(phi)*np.cos(d_)*np.cos(w)
    tz = np.degrees(np.arccos(np.clip(cz, -1, 1)))
    if tz > 80: continue
    sz = np.sin(np.radians(tz))
    cg = np.clip((cz*np.sin(phi) - np.sin(d_))/(sz*np.cos(phi)), -1, 1)
    az = np.degrees(np.sign(ha)*abs(np.arccos(cg))) + 180.0
    ref = float(np.clip(np.cos(np.radians(pvlib.irradiance.aoi(
        surface_tilt=beta, surface_azimuth=0.0, solar_zenith=tz, solar_azimuth=az))), 0, 1))
    worst_pv = max(worst_pv, abs(cos_aoi(lat, beta, dec, ha) - ref)); n += 1
check(f"ecart max vs pvlib = {worst_pv:.3e} sur {n} points australs", worst_pv < 1e-9)

print(); print("="*78); print("CRITERE 3 -- latitude nulle"); print("="*78)
b = 20.0
at0, flat = cos_aoi(0.0, b, 0.0, 0.0), cos_aoi(0.0, 0.0, 0.0, 0.0)
check(f"cos(theta) a lat=0 = {at0:.6f}, attendu cos(20) = {np.cos(np.radians(b)):.6f}",
      abs(at0 - np.cos(np.radians(b))) < 1e-12)
check("l'inclinaison n'est pas annulee a lat=0", abs(at0 - flat) > 1e-6)

print(); print("="*78); print("CRITERE 4 -- zero contamination (memes coordonnees)"); print("="*78)
df = pd.read_csv(REPO/"outputs"/"csv"/"hourly_Ouagadougou.csv")
df["time_UTC"] = pd.to_datetime(df["time_UTC"]); df = df.set_index("time_UTC")
X = [16.37, 3.026, 1.501, 2.64, 31.8, 0.0502]; KEY = "_runtime_12.37_-1.53"
def site(p, y):
    s = copy.deepcopy(dict(config_00.SITES["Ouagadougou"])); s["PAR_sat"]=p; s["fruit_yield_t_ha"]=y; return s
cfgs = {"tomate": (174.0, 60.0), "mais": (400.0, 12.0)}
serial = {k: evaluate_point(KEY, site(*v), df, X, scenario="S0", PAR_sat=v[0]) for k, v in cfgs.items()}
res, errs = {}, {}
jobs = [(f"{k}#{i}", site(*v), v[0]) for i in range(12) for k, v in cfgs.items()]
bar = threading.Barrier(len(jobs))
def w(label, s_, p_):
    try:
        bar.wait(); res[label] = evaluate_point(KEY, s_, df, X, scenario="S0", PAR_sat=p_)
    except Exception as e: errs[label] = f"{type(e).__name__}: {e}"
ts = [threading.Thread(target=w, args=j) for j in jobs]
[t.start() for t in ts]; [t.join() for t in ts]
bad = [l for l, r in res.items() if abs(r["eLER_defA"] - serial[l.split("#")[0]]["eLER_defA"]) > 1e-12]
check(f"{len(jobs)} threads simultanes, {len(errs)} erreurs, {len(bad)} contamines",
      not errs and not bad, f"{errs if errs else ''}{bad if bad else ''}")

print(); print("="*78); print("CRITERE 5 -- aucune dependance MATLAB"); print("="*78)
mods = set()
for f in list((REPO/"notebooks").glob("*.py")) + list((REPO/"webapp"/"backend").glob("*.py")) + list((REPO/"tests").glob("*.py")):
    t = f.read_text(encoding="utf-8", errors="replace").lower()
    if "matlab" in t and "matlab/" not in t: mods.add(f.name)
check("aucun module Python n'importe MATLAB", not mods, str(mods))

print(); print("="*78); print("CRITERE 6 -- erreurs de configuration explicites"); print("="*78)
try:
    _ = config_00.SITES["_inconnu"]; check("UnknownSiteError levee", False)
except UnknownSiteError as e:
    check("UnknownSiteError, message lisible", "Unknown site" in str(e))
    check("compatible avec except KeyError", issubclass(UnknownSiteError, KeyError))

print(); print("="*78); print("CRITERE 7 -- notebook et export en phase (drift CI)"); print("="*78)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from regen import regenerate
got = regenerate(Path("/tmp/regen_check.py"))
ref_txt = (REPO/"notebooks"/"simulation_functions.py").read_text(encoding="utf-8")
check("simulation_functions.py == export du notebook", got == ref_txt)

print(); print("="*78)
print("RESULTAT :", "TOUS LES CRITERES PASSENT" if ok_all else "ECHEC")
print("="*78)
sys.exit(0 if ok_all else 1)
