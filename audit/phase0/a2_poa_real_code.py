"""
AUDIT PHASE 0 -- A2/A3 : impact sur le VRAI code, et bornes d'inclinaison
========================================================================

LECTURE SEULE. Aucun fichier du depot n'est ecrit. La correction est
appliquee par monkeypatch EN MEMOIRE, uniquement pour mesurer l'ecart.

Protocole
---------
PVGIS n'est pas joignable depuis ce conteneur, ce qui impose -- et permet --
un protocole plus propre : on garde UN SEUL climat reel (la TMY de Konya,
8760 h) et on fait varier la seule latitude. Toute difference observee est
donc imputable a la geometrie et a rien d'autre. Un jeu de donnees austral
reel melangerait climat et geometrie ; ici les deux sont separes.

Limite assumee : le climat de Konya n'est pas celui de Lusaka. Les valeurs
absolues de POA ne sont pas representatives ; seuls les ECARTS RELATIFS
entre code actuel et code corrige, a latitude donnee, le sont.

A2 : ecart de POA et de productible PV, par latitude, sur le vrai pipeline.
A3 : optimum radiatif d'inclinaison vs les bornes reellement imposees
     (DE_SETTINGS + min_beta_deg), et verification sur les 4 sites publies.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent / "repo"
sys.path.insert(0, str(REPO / "notebooks"))

import config_00                      # noqa: E402  (non modifie)
import simulation_functions as sim    # noqa: E402  (non modifie)

OUT = "a2_poa_real_code"
P = config_00.PARAMS


# ----------------------------------------------------------------------
# Version corrigee de compute_POA : identique a l'originale SAUF phi_eff.
# Appliquee par monkeypatch, jamais ecrite sur disque.
# ----------------------------------------------------------------------
def compute_POA_corrected(df, beta_deg, lat_deg):
    beta_r = np.radians(beta_deg)
    lat_r = np.radians(lat_deg)

    dec_r = np.radians(df["declination_deg"].values)
    ha_r = np.radians(df["hour_angle_deg"].values)

    # SEULE DIFFERENCE : le capteur se tourne vers l'equateur.
    #   phi - beta  dans l'hemisphere nord,  phi + beta  dans l'hemisphere sud
    phi_eff = lat_r - np.sign(lat_r) * beta_r

    cos_theta_i = np.clip(
        np.sin(dec_r) * np.sin(phi_eff)
        + np.cos(dec_r) * np.cos(ha_r) * np.cos(phi_eff),
        0, 1,
    )

    GHI = df["GHI"].values
    DNI = df["DNI"].values
    DHI = df["DHI"].values

    G_beam = DNI * cos_theta_i
    G_diffuse = DHI * (1 + np.cos(beta_r)) / 2
    G_reflect = GHI * P["rho_ground"] * (1 - np.cos(beta_r)) / 2
    G_front = np.clip(G_beam + G_diffuse + G_reflect, 0, None)

    phi_rear = 0.12
    G_rear = np.clip(GHI * P["rho_ground"] * (1 - df["F_shad"].values) * phi_rear, 0, None)

    df = df.copy()
    df["G_front"] = G_front
    df["G_rear"] = G_rear
    df["G_eff"] = G_front + P["bifaciality"] * G_rear
    return df


def load_konya():
    df = pd.read_csv(REPO / "outputs" / "csv" / "hourly_Konya.csv")
    df["time_UTC"] = pd.to_datetime(df["time_UTC"])
    return df.set_index("time_UTC")


def annual_poa(df_raw, lat, lon, beta, d_row, H_m, poa_fn):
    """POA front et productible PV annuels pour une latitude donnee."""
    d = sim.compute_solar_angles(df_raw, lat, lon)
    d = sim.compute_shading(d, beta, d_row, H_m)
    d = poa_fn(d, beta, lat)
    n_mod = sim.modules_per_ha(beta, d_row)
    d = sim.compute_PV(d, n_mod, 0.0)
    return (float(d["G_front"].sum()) / 1000,          # kWh/m2/an
            float(d["E_PV_kWh"].sum()) / 1000,         # MWh/ha/an
            float(d["F_shad"].mean()),
            float(d["solar_elevation_deg"].max()))


# ======================================================================
# A2 : ecart par latitude sur le vrai pipeline
# ======================================================================
def a2_latitude_transect(df_raw):
    lats = np.arange(-35, 38, 2.5)
    beta, d_row, H_m = 20.0, 6.0, 2.5       # design fixe : seule la latitude varie
    rows = []
    for lat in lats:
        poa_c, pv_c, fsh_c, elev = annual_poa(df_raw, lat, 30.0, beta, d_row, H_m,
                                               sim.compute_POA)
        poa_f, pv_f, fsh_f, _ = annual_poa(df_raw, lat, 30.0, beta, d_row, H_m,
                                            compute_POA_corrected)
        rows.append({
            "lat": lat,
            "POA_actuel_kWh_m2": round(poa_c, 1),
            "POA_corrige_kWh_m2": round(poa_f, 1),
            "ecart_POA_pct": round(100 * (poa_c - poa_f) / poa_f, 2),
            "PV_actuel_MWh_ha": round(pv_c, 1),
            "PV_corrige_MWh_ha": round(pv_f, 1),
            "ecart_PV_pct": round(100 * (pv_c - pv_f) / pv_f, 2),
            "F_shad_identique": abs(fsh_c - fsh_f) < 1e-12,
        })
    return pd.DataFrame(rows)


# ======================================================================
# A3 : optimum radiatif d'inclinaison vs bornes imposees
# ======================================================================
def a3_tilt_optimum(df_raw):
    lats = np.arange(-35, 38, 2.5)
    betas = np.arange(0, 46, 1.0)
    d_row, H_m = 6.0, 2.5
    lo_de = config_00.DE_SETTINGS["bounds"][0][0]      # borne basse beta du DE

    rows = []
    for lat in lats:
        best_b, best_poa = None, -1.0
        for b in betas:
            poa, _, _, _ = annual_poa(df_raw, lat, 30.0, b, d_row, H_m,
                                       compute_POA_corrected)
            if poa > best_poa:
                best_poa, best_b = poa, b
        # POA atteignable si la borne basse du DE s'impose
        poa_at_bound, _, _, _ = annual_poa(df_raw, lat, 30.0, max(best_b, lo_de),
                                            d_row, H_m, compute_POA_corrected)
        rows.append({
            "lat": lat,
            "beta_opt_radiatif": best_b,
            "regle_manuscrit_0.45|phi|": round(0.45 * abs(lat), 1),
            "borne_basse_DE": lo_de,
            "borne_active": best_b < lo_de,
            "perte_si_borne_pct": round(100 * (poa_at_bound - best_poa) / best_poa, 2),
        })
    return pd.DataFrame(rows)


def a3_published_sites():
    """Les 4 sites publies convergent-ils sur leur borne basse d'inclinaison ?"""
    import json
    opt = json.loads((REPO / "outputs" / "csv" / "optimal_params.json").read_text())
    lo_de = config_00.DE_SETTINGS["bounds"][0][0]
    rows = []
    for name, cfg in config_00.SITES.items():
        if name not in opt:
            continue
        b_opt = opt[name]["beta_deg"]
        floor = max(lo_de, cfg["min_beta_deg"])
        rows.append({
            "site": name,
            "lat": cfg["lat"],
            "beta_optimise": round(b_opt, 2),
            "min_beta_deg": cfg["min_beta_deg"],
            "borne_effective": floor,
            "marge_au_dessus_borne": round(b_opt - floor, 2),
            "regle_0.45|phi|": round(0.45 * abs(cfg["lat"]), 1),
        })
    return pd.DataFrame(rows)


# ======================================================================
if __name__ == "__main__":
    pd.set_option("display.width", 150)
    pd.set_option("display.max_rows", 200)

    df_raw = load_konya()
    print(f"Climat de reference : Konya TMY, {len(df_raw)} heures, "
          f"GHI {df_raw['GHI'].sum()/1000:.0f} kWh/m2/an")
    print("Seule la latitude varie -- le climat est constant, donc tout ecart")
    print("est imputable a la geometrie d'orientation.\n")

    print("=" * 96)
    print("A2  Ecart de POA et de productible PV par latitude (vrai pipeline)")
    print("=" * 96)
    t = a2_latitude_transect(df_raw)
    print(t.to_string(index=False))
    t.to_csv(f"{OUT}_transect.csv", index=False)

    south = t[t.lat < 0]
    print()
    print(f"  Hemisphere NORD : ecart POA max {t[t.lat > 0].ecart_POA_pct.abs().max():.3f} %")
    print(f"  Hemisphere SUD  : ecart POA de {south.ecart_POA_pct.min():.1f} % "
          f"a {south.ecart_POA_pct.max():.1f} %")
    print(f"  F_shad identique partout : {bool(t.F_shad_identique.all())} "
          f"(confirme que compute_shading ne depend pas de l'azimut)")

    print()
    print("=" * 96)
    print("A3  Optimum radiatif d'inclinaison vs bornes imposees")
    print("=" * 96)
    o = a3_tilt_optimum(df_raw)
    print(o.to_string(index=False))
    o.to_csv(f"{OUT}_tilt.csv", index=False)

    binding = o[o.borne_active]
    print()
    print(f"  Latitudes ou la borne basse du DE ({config_00.DE_SETTINGS['bounds'][0][0]} deg) "
          f"est active : {len(binding)} / {len(o)}")
    if len(binding):
        print(f"  Plage concernee : {binding.lat.min():.1f} deg a {binding.lat.max():.1f} deg")
        print(f"  Perte de POA maximale imputable a la borne : "
              f"{binding.perte_si_borne_pct.min():.2f} %")

    print()
    print("=" * 96)
    print("A3b  Les 4 sites publies convergent-ils sur leur borne ?")
    print("=" * 96)
    s = a3_published_sites()
    print(s.to_string(index=False))
    s.to_csv(f"{OUT}_sites.csv", index=False)
