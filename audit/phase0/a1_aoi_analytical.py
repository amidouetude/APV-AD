"""
AUDIT PHASE 0 -- A1 : angle d'incidence, comparaison analytique
==============================================================

LECTURE SEULE. Ne modifie aucun fichier du depot.

Objet
-----
Isoler l'erreur d'ORIENTATION (hypothese plein sud) de toute autre source
d'ecart, en comparant trois calculs de cos(theta_i) sur exactement la meme
position solaire analytique :

  1. FORMULE ACTUELLE  -- celle de compute_POA :
         cos t = sin d . sin(phi - beta) + cos d . cos w . cos(phi - beta)
     C'est la forme fermee de Duffie & Beckman pour gamma = 0 (plein sud).

  2. FORMULE CORRIGEE  -- Duffie & Beckman complete, capteur tourne vers
     l'equateur : plein sud si phi > 0, plein nord si phi < 0.

  3. PVLIB            -- oracle independant, alimente avec la MEME position
     solaire analytique, donc tout ecart entre (2) et (3) est une erreur
     d'implementation de ma part, pas un ecart de modele solaire.

Conventions d'azimut (piege classique) :
  - Duffie & Beckman : gamma = 0 au SUD, positif vers l'ouest.
  - pvlib            : azimut = 0 au NORD, sens horaire.
  => azimut_pvlib = gamma_DB + 180
"""
import numpy as np
import pandas as pd
import pvlib

OUT = "a1_aoi_analytical"


# ----------------------------------------------------------------------
# Position solaire analytique (independante du depot, pour isoler la geometrie)
# ----------------------------------------------------------------------
def solar_position(phi_deg, dec_deg, omega_deg):
    """Zenith solaire et azimut (convention pvlib : 0 = Nord, horaire)."""
    phi = np.radians(phi_deg)
    dec = np.radians(dec_deg)
    w = np.radians(omega_deg)

    cos_z = np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.cos(w)
    cos_z = np.clip(cos_z, -1, 1)
    theta_z = np.arccos(cos_z)

    sin_z = np.sin(theta_z)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_gs = (cos_z * np.sin(phi) - np.sin(dec)) / (sin_z * np.cos(phi))
    cos_gs = np.clip(cos_gs, -1, 1)
    gamma_s = np.sign(omega_deg) * np.abs(np.arccos(cos_gs))   # depuis le sud, + ouest
    az_pvlib = np.degrees(gamma_s) + 180.0
    return np.degrees(theta_z), az_pvlib


# ----------------------------------------------------------------------
# Les trois calculs d'angle d'incidence
# ----------------------------------------------------------------------
def aoi_current_code(phi_deg, beta_deg, dec_deg, omega_deg):
    """Exactement la ligne de compute_POA (plein sud implicite)."""
    phi = np.radians(phi_deg)
    beta = np.radians(beta_deg)
    dec = np.radians(dec_deg)
    w = np.radians(omega_deg)
    return np.clip(
        np.sin(dec) * np.sin(phi - beta)
        + np.cos(dec) * np.cos(w) * np.cos(phi - beta),
        0, 1,
    )


def aoi_duffie_full(phi_deg, beta_deg, gamma_deg, dec_deg, omega_deg):
    """Duffie & Beckman eq. 1.6.5, azimut de surface gamma explicite (0 = sud)."""
    phi = np.radians(phi_deg)
    beta = np.radians(beta_deg)
    g = np.radians(gamma_deg)
    dec = np.radians(dec_deg)
    w = np.radians(omega_deg)
    cos_t = (
        np.sin(dec) * np.sin(phi) * np.cos(beta)
        - np.sin(dec) * np.cos(phi) * np.sin(beta) * np.cos(g)
        + np.cos(dec) * np.cos(phi) * np.cos(beta) * np.cos(w)
        + np.cos(dec) * np.sin(phi) * np.sin(beta) * np.cos(g) * np.cos(w)
        + np.cos(dec) * np.sin(beta) * np.sin(g) * np.sin(w)
    )
    return np.clip(cos_t, 0, 1)


def equator_facing_gamma(phi_deg):
    """gamma (convention D&B) d'un capteur tourne vers l'equateur."""
    return np.where(np.asarray(phi_deg) >= 0, 0.0, 180.0)


# ======================================================================
# 1. Verification croisee de MON implementation contre pvlib
# ======================================================================
def check_against_pvlib():
    rows = []
    rng = np.random.default_rng(42)
    for _ in range(4000):
        phi = rng.uniform(-40, 40)
        beta = rng.uniform(0, 45)
        gamma = rng.uniform(-180, 180)
        dec = rng.uniform(-23.45, 23.45)
        omega = rng.uniform(-90, 90)

        theta_z, az_pv = solar_position(phi, dec, omega)
        if theta_z > 85:            # trop bas : sin(z) ~ 0, azimut mal conditionne
            continue

        mine = aoi_duffie_full(phi, beta, gamma, dec, omega)
        ref = np.cos(np.radians(
            pvlib.irradiance.aoi(
                surface_tilt=beta,
                surface_azimuth=gamma + 180.0,
                solar_zenith=theta_z,
                solar_azimuth=az_pv,
            )
        ))
        ref = np.clip(ref, 0, 1)
        rows.append((phi, beta, gamma, dec, omega, mine, ref, abs(mine - ref)))

    df = pd.DataFrame(rows, columns=["phi", "beta", "gamma", "dec", "omega",
                                      "mine", "pvlib", "abs_err"])
    return df


# ======================================================================
# 2. Erreur d'orientation par latitude, au midi solaire
# ======================================================================
def noon_error_by_latitude():
    lats = np.arange(-35, 38, 2.5)
    rows = []
    for phi in lats:
        beta = 0.45 * abs(phi)          # Finding 2 du manuscrit : 0.4-0.5 x |phi|
        gamma = float(equator_facing_gamma(phi))
        for label, dec in [("equinoxe", 0.0),
                            ("solstice juin", 23.45),
                            ("solstice decembre", -23.45)]:
            cur = float(aoi_current_code(phi, beta, dec, 0.0))
            cor = float(aoi_duffie_full(phi, beta, gamma, dec, 0.0))
            rows.append({
                "lat": phi, "beta": round(beta, 2), "saison": label,
                "cos_ti_actuel": round(cur, 4),
                "cos_ti_corrige": round(cor, 4),
                "ecart_rel_pct": round(100 * (cur - cor) / cor, 2) if cor > 1e-9 else np.nan,
            })
    return pd.DataFrame(rows)


# ======================================================================
# 3. Erreur integree sur une journee claire (proxy de l'impact energetique)
# ======================================================================
def daily_integrated_error():
    """
    Somme de cos(theta_i) sur les heures de jour, pour 12 declinaisons
    reparties sur l'annee. Proxy geometrique pur du rayonnement direct
    recu : aucune donnee meteo, donc aucune autre source d'ecart.
    """
    lats = np.arange(-35, 38, 2.5)
    decs = [np.degrees(0.4093 * np.sin(2 * np.pi * (284 + d) / 365))
            for d in (15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)]
    omegas = np.arange(-180, 180, 1.0)      # pas de 4 min

    rows = []
    for phi in lats:
        beta = 0.45 * abs(phi)
        gamma = float(equator_facing_gamma(phi))
        tot_cur = tot_cor = 0.0
        for dec in decs:
            theta_z, _ = solar_position(phi, dec, omegas)
            day = theta_z < 90
            tot_cur += float(np.sum(aoi_current_code(phi, beta, dec, omegas)[day]))
            tot_cor += float(np.sum(aoi_duffie_full(phi, beta, gamma, dec, omegas)[day]))
        rows.append({
            "lat": phi, "beta": round(beta, 2),
            "somme_cos_ti_actuel": round(tot_cur, 1),
            "somme_cos_ti_corrige": round(tot_cor, 1),
            "ecart_rel_pct": round(100 * (tot_cur - tot_cor) / tot_cor, 2),
        })
    return pd.DataFrame(rows)


# ======================================================================
if __name__ == "__main__":
    pd.set_option("display.width", 130)
    pd.set_option("display.max_rows", 200)

    print("=" * 74)
    print("A1.1  Verification de mon implementation Duffie & Beckman vs pvlib")
    print("=" * 74)
    chk = check_against_pvlib()
    print(f"  echantillons retenus       : {len(chk)}")
    print(f"  ecart absolu max           : {chk.abs_err.max():.3e}")
    print(f"  ecart absolu median        : {chk.abs_err.median():.3e}")
    ok = chk.abs_err.max() < 1e-9
    print(f"  -> {'CONFORME' if ok else 'ECART -- implementation a revoir'}")
    chk.to_csv(f"{OUT}_pvlib_crosscheck.csv", index=False)

    print()
    print("=" * 74)
    print("A1.2  Erreur d'orientation au midi solaire, par latitude")
    print("=" * 74)
    noon = noon_error_by_latitude()
    eq = noon[noon.saison == "equinoxe"]
    print(eq.to_string(index=False))
    noon.to_csv(f"{OUT}_noon.csv", index=False)

    print()
    print("=" * 74)
    print("A1.3  Erreur integree sur l'annee (geometrie pure, sans meteo)")
    print("=" * 74)
    daily = daily_integrated_error()
    print(daily.to_string(index=False))
    daily.to_csv(f"{OUT}_daily.csv", index=False)

    print()
    print("-" * 74)
    south = daily[daily.lat < 0]
    north = daily[daily.lat > 0]
    print(f"Hemisphere NORD  : ecart max {north.ecart_rel_pct.abs().max():.2f} %")
    print(f"Hemisphere SUD   : ecart max {south.ecart_rel_pct.abs().max():.2f} %, "
          f"median {south.ecart_rel_pct.median():.2f} %")
    print("-" * 74)
