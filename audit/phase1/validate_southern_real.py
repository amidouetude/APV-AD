"""
PHASE 1 -- validation australe sur donnees PVGIS reelles
========================================================

A executer DEPUIS LE POSTE UTILISATEUR : le conteneur d'audit n'a pas
d'acces reseau a PVGIS (proxy, HTTP 403), ce qui avait impose en phase 0 un
protocole a climat constant. Ce script fait ce que ce protocole ne pouvait
pas faire : verifier le comportement INTEGRE sur des chroniques climatiques
australes reelles, avec leur saisonnalite propre.

Ce qui est verifie
------------------
  1. Le site de controle NORD reproduit exactement les valeurs publiees.
  2. POA et productible PV, avant / apres correction, sur climat reel.
  3. L'orientation retenue est bien tournee vers l'equateur.
  4. L'optimum radiatif d'inclinaison en climat austral reel -- la grandeur
     que la phase 0 avait explicitement marquee "non transferable", faute
     de saisonnalite australe.
  5. Un run complet de l'optimiseur DE en hemisphere sud (option --de).

Les TMY telechargees et toutes les sorties sont conservees comme artefacts
d'audit dans audit/phase1/southern_validation/.

Usage
-----
    python validate_southern_real.py --probe     # test d'acces PVGIS seul
    python validate_southern_real.py             # validation complete
    python validate_southern_real.py --de        # + un run DE complet (~200 s)
"""
import argparse
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]                      # .../APV_AD_Project
OUT = HERE / "southern_validation"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO / "notebooks"))
sys.path.insert(0, str(REPO / "webapp" / "backend"))

LOG = OUT / "run.log"


def log(msg=""):
    line = str(msg)
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as f:
        f.write(line + "\n")


# ----------------------------------------------------------------------
# Sites : un controle nord, puis un transect austral
# ----------------------------------------------------------------------
SITES = [
    # nom,           lat,     lon,    saison (MM-DD -> MM-DD), role
    ("Ouagadougou",  12.37,  -1.53,  ("07-01", "10-30"), "controle NORD"),
    ("Nairobi",      -1.29,  36.82,  ("03-15", "07-15"), "quasi-equatorial"),
    ("Lusaka",      -15.42,  28.28,  ("11-15", "04-15"), "austral tropical"),
    ("Windhoek",    -22.57,  17.08,  ("11-01", "04-30"), "austral subtropical"),
    ("Le Cap",      -33.92,  18.42,  ("09-01", "12-31"), "austral tempere"),
]

# Les dates de saison ci-dessus sont des HYPOTHESES plausibles (semis en
# debut de saison des pluies dans l'hemisphere sud), PAS un calendrier
# cultural source FAO/GIEWS. Elles n'influencent ni POA ni le productible
# PV -- seulement les volets culture et eau. Marquees comme telles partout.
SEASON_SOURCE = "hypothese -- PAS un calendrier cultural source (FAO/GIEWS a faire)"

X_DEFAULT = [16.37, 3.026, 1.501, 2.64, 31.8, 0.0502]   # design d'Ouagadougou


# ----------------------------------------------------------------------
def legacy_cos_aoi(dec_deg, ha_deg, lat_deg, beta_deg):
    """
    Forme PLEIN SUD d'avant la correction, reimplantee ici pour le
    comparatif avant/apres. Rien n'est modifie dans le depot.
    """
    dec_r = np.radians(dec_deg)
    ha_r = np.radians(ha_deg)
    phi_minus_beta = np.radians(lat_deg - beta_deg)
    return np.clip(
        np.sin(dec_r) * np.sin(phi_minus_beta)
        + np.cos(dec_r) * np.cos(ha_r) * np.cos(phi_minus_beta),
        0, 1,
    )


def poa_with(df, beta, lat, legacy=False):
    """POA front annuel (kWh/m2) et productible PV (MWh/ha), corrige ou legacy."""
    import config_00
    import simulation_functions as sim

    P = config_00.PARAMS
    d = sim.compute_solar_angles(df, lat, df.attrs["lon"])
    d = sim.compute_shading(d, beta, 6.0, 2.5)

    if legacy:
        beta_r = np.radians(beta)
        cos_ti = legacy_cos_aoi(d["declination_deg"].values,
                                 d["hour_angle_deg"].values, lat, beta)
        G_front = np.clip(
            d["DNI"].values * cos_ti
            + d["DHI"].values * (1 + np.cos(beta_r)) / 2
            + d["GHI"].values * P["rho_ground"] * (1 - np.cos(beta_r)) / 2,
            0, None)
        G_rear = np.clip(d["GHI"].values * P["rho_ground"]
                          * (1 - d["F_shad"].values) * 0.12, 0, None)
        d = d.copy()
        d["G_front"] = G_front
        d["G_rear"] = G_rear
        d["G_eff"] = G_front + P["bifaciality"] * G_rear
    else:
        d = sim.compute_POA(d, beta, lat)

    n_mod = sim.modules_per_ha(beta, 6.0)
    d = sim.compute_PV(d, n_mod, 0.0)
    return float(d["G_front"].sum()) / 1000, float(d["E_PV_kWh"].sum()) / 1000


# ----------------------------------------------------------------------
def fetch(lat, lon, name, season):
    from pvgis_client import fetch_tmy, add_growing_season
    cache = OUT / f"tmy_{name.replace(' ', '_')}.csv"
    if cache.exists():
        log(f"    TMY en cache : {cache.name}")
        df = pd.read_csv(cache)
        df["time_UTC"] = pd.to_datetime(df["time_UTC"])
        df = df.set_index("time_UTC")
    else:
        t0 = time.time()
        df = fetch_tmy(lat, lon)
        log(f"    PVGIS {time.time() - t0:.1f} s | db={df.attrs['radiation_db']} "
            f"| alt={df.attrs['elevation_m']} m")
        df.to_csv(cache)
        (OUT / f"tmy_{name.replace(' ', '_')}.meta.json").write_text(json.dumps({
            "lat": lat, "lon": lon, "radiation_db": df.attrs["radiation_db"],
            "elevation_m": df.attrs["elevation_m"],
            "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }, indent=2))
    df = add_growing_season(df, season[0], season[1], source=SEASON_SOURCE)
    df.attrs["lat"], df.attrs["lon"] = lat, lon
    return df


# ======================================================================
def main(run_de=False):
    import config_00

    log("=" * 86)
    log("VALIDATION AUSTRALE -- donnees PVGIS reelles")
    log(f"depot : {REPO}")
    log("=" * 86)

    rows, sweeps = [], []
    for name, lat, lon, season, role in SITES:
        log(f"\n--- {name} ({lat:+.2f}, {lon:+.2f})  [{role}] ---")
        try:
            df = fetch(lat, lon, name, season)
        except Exception as exc:
            log(f"    ECHEC PVGIS : {type(exc).__name__}: {exc}")
            rows.append({"site": name, "lat": lat, "erreur": str(exc)})
            continue

        ghi = df["GHI"].sum() / 1000
        beta = 20.0
        poa_c, pv_c = poa_with(df, beta, lat, legacy=False)
        poa_l, pv_l = poa_with(df, beta, lat, legacy=True)

        orient = "nord (equateur)" if lat < 0 else "sud (equateur)"
        rows.append({
            "site": name, "role": role, "lat": lat, "lon": lon,
            "GHI_kWh_m2": round(ghi, 1),
            "orientation": orient,
            "POA_corrige": round(poa_c, 1), "POA_legacy": round(poa_l, 1),
            "ecart_POA_pct": round(100 * (poa_l - poa_c) / poa_c, 2),
            "PV_corrige_MWh_ha": round(pv_c, 1), "PV_legacy_MWh_ha": round(pv_l, 1),
            "ecart_PV_pct": round(100 * (pv_l - pv_c) / pv_c, 2),
        })
        log(f"    GHI {ghi:.0f} kWh/m2 | orientation {orient}")
        log(f"    POA corrige {poa_c:7.1f}  legacy {poa_l:7.1f}  "
            f"ecart {100*(poa_l-poa_c)/poa_c:+6.2f} %")
        log(f"    PV  corrige {pv_c:7.1f}  legacy {pv_l:7.1f}  "
            f"ecart {100*(pv_l-pv_c)/pv_c:+6.2f} %")

        # Optimum radiatif d'inclinaison, CLIMAT REEL
        best_b, best_poa = None, -1.0
        for b in np.arange(0, 46, 1.0):
            p, _ = poa_with(df, float(b), lat)
            if p > best_poa:
                best_poa, best_b = p, float(b)
        floor = max(config_00.DE_SETTINGS["bounds"][0][0], 15.0)
        sweeps.append({
            "site": name, "lat": lat,
            "beta_opt_radiatif": best_b,
            "regle_0.45|phi|": round(0.45 * abs(lat), 1),
            "borne_basse_DE": floor,
            "borne_active": best_b < floor,
            "POA_a_l_optimum": round(best_poa, 1),
            "POA_a_la_borne": round(poa_with(df, max(best_b, floor), lat)[0], 1),
        })
        log(f"    beta optimal (climat reel) = {best_b:.0f} deg "
            f"| regle 0.45|phi| = {0.45*abs(lat):.1f} | borne DE = {floor:.0f}"
            f" | borne active : {best_b < floor}")

    df_rows = pd.DataFrame(rows)
    df_sweeps = pd.DataFrame(sweeps)
    df_rows.to_csv(OUT / "results_poa_pv.csv", index=False)
    df_sweeps.to_csv(OUT / "results_tilt_optimum.csv", index=False)

    log("\n" + "=" * 86)
    log("SYNTHESE -- POA / PV")
    log("=" * 86)
    log(df_rows.to_string(index=False))
    log("\n" + "=" * 86)
    log("SYNTHESE -- optimum d'inclinaison en climat reel")
    log("=" * 86)
    log(df_sweeps.to_string(index=False))

    # -- Controle nord : ecart doit etre nul --
    north = df_rows[df_rows.lat > 0]
    if len(north):
        worst = north["ecart_POA_pct"].abs().max()
        log(f"\n[{'OK ' if worst < 1e-9 else 'FAIL'}] controle NORD : "
            f"ecart legacy/corrige = {worst:.2e} % (doit etre nul)")

    south = df_rows[df_rows.lat < 0]
    if len(south):
        log(f"[INFO] hemisphere SUD : le code d'avant sous-estimait le POA de "
            f"{south['ecart_POA_pct'].min():.1f} % a {south['ecart_POA_pct'].max():.1f} %")

    # ------------------------------------------------------------------
    if run_de:
        log("\n" + "=" * 86)
        log("RUN DE COMPLET -- hemisphere sud (Lusaka)")
        log("=" * 86)
        from optimizer import optimize_site
        import copy
        name, lat, lon, season, _ = SITES[2]          # Lusaka
        df = fetch(lat, lon, name, season)
        site = copy.deepcopy(dict(config_00.SITES["Ouagadougou"]))
        site.update({"lat": lat, "lon": lon,
                     "elevation_m": 1280, "min_beta_deg": 12.0})
        t0 = time.time()
        gens = []
        x_opt, full, hist = optimize_site(
            f"_val_{name}", site, df, scenario="S0", PAR_sat=174.0,
            progress_cb=lambda g, b: gens.append((g, b)))
        log(f"    duree {time.time() - t0:.0f} s, {len(hist)} generations")
        log(f"    beta={x_opt[0]:.2f}  d_row={x_opt[1]:.2f}  H={x_opt[2]:.2f}")
        log(f"    V_dig={x_opt[3]:.2f}  HRT={x_opt[4]:.1f}  f_PV_heat={x_opt[5]:.4f}")
        log(f"    eLER = {full['eLER_defA']:.4f} | LER_crop {full['LER_crop']:.3f} "
            f"| LER_PV {full['LER_PV_defA']:.3f} | LER_biogas {full['LER_biogas']:.3f}")
        log(f"    PV {full['PV_total_MWh_ha']:.1f} MWh/ha | "
            f"biogaz {full['biogas_total_MWh_ha']:.2f} MWh/ha")
        (OUT / "de_lusaka.json").write_text(json.dumps({
            "site": name, "lat": lat, "lon": lon,
            "x_opt": list(map(float, x_opt)),
            "result": {k: v for k, v in full.items() if isinstance(v, (int, float, str))},
            "history": [float(h) for h in hist],
            "season_source": SEASON_SOURCE,
        }, indent=2, default=str))
        log(f"    artefact : {(OUT / 'de_lusaka.json').name}")

    log("\n" + "=" * 86)
    log(f"TERMINE -- artefacts dans {OUT}")
    log("=" * 86)


# ======================================================================
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true", help="tester l'acces PVGIS seul")
    ap.add_argument("--de", action="store_true", help="ajouter un run DE complet")
    a = ap.parse_args()

    LOG.write_text("", encoding="utf-8")
    try:
        if a.probe:
            from pvgis_client import fetch_tmy
            t0 = time.time()
            d = fetch_tmy(-15.42, 28.28)
            log(f"PROBE OK -- PVGIS joignable, {len(d)} lignes en "
                f"{time.time()-t0:.1f} s, GHI {d['GHI'].sum()/1000:.0f} kWh/m2")
            log(f"python {sys.version.split()[0]} | pandas {pd.__version__} "
                f"| numpy {np.__version__}")
        else:
            main(run_de=a.de)
    except Exception:
        log("ECHEC :\n" + traceback.format_exc())
        sys.exit(1)
