"""
AUDIT PHASE 0 -- C : gel des references, et D : carte des hypotheses d'orientation
=================================================================================

LECTURE SEULE sur le depot. Ecrit uniquement ses propres fichiers de sortie
dans le dossier d'audit.

C : rejoue le chemin de code ACTUEL pour les 4 sites publies et fige toutes
    les sorties dans un JSON. C'est la reference de non-regression : apres
    correction, ces valeurs doivent etre identiques (hemisphere nord => la
    correction ne doit rien changer du tout).

D : recense chaque endroit du depot ou une orientation plein sud, un azimut
    ou une borne d'inclinaison est supposee implicitement.
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent / "repo"
sys.path.insert(0, str(REPO / "notebooks"))
sys.path.insert(0, str(REPO / "webapp" / "backend"))

import config_00                      # noqa: E402
from model_bridge import evaluate_catalog_site  # noqa: E402
from economics import compute_4e      # noqa: E402

TOL_REL = 1e-9          # tolerance de non-regression : identite stricte attendue


# ======================================================================
# C : gel des references
# ======================================================================
def freeze_references():
    opt = json.loads((REPO / "outputs" / "csv" / "optimal_params.json").read_text())
    published = pd.read_csv(REPO / "outputs" / "csv" / "results_summary.csv")

    frozen, compare = {}, []
    for site in config_00.SITES:
        if site not in opt:
            continue
        o = opt[site]
        x = [o["beta_deg"], o["d_row_m"], o["H_m_m"],
             o["V_dig_m3"], o["HRT_days"], o["f_PV_heat"]]

        res = evaluate_catalog_site(site, x, scenario="S0")
        econ = compute_4e(config_00.SITES[site], res, scenario="S0")

        frozen[site] = {
            "design_vars": x,
            "scenario": "S0",
            "result": {k: v for k, v in res.items()
                       if isinstance(v, (int, float, str))},
            "economics": econ,
        }

        row = published[(published["site"] == site) & (published["scenario"] == "S0")]
        for field in ("eLER_defA", "LER_crop", "LER_PV_defA", "LER_biogas",
                      "PV_total_MWh_ha", "biogas_total_MWh_ha",
                      "W_saved_mm_season", "LER_water_season_pct"):
            if field in res and len(row) and field in row.columns:
                got, exp = float(res[field]), float(row.iloc[0][field])
                compare.append({
                    "site": site, "champ": field,
                    "rejoue": round(got, 4), "publie": round(exp, 4),
                    "ecart_pct": round(100 * (got - exp) / exp, 4) if exp else None,
                })
    return frozen, pd.DataFrame(compare)


# ======================================================================
# D : carte des hypotheses d'orientation
# ======================================================================
PATTERNS = [
    (r"lat_r\s*-\s*beta_r",            "forme plein sud codee en dur (Duffie gamma=0)"),
    (r"delta_gamma\s*=\s*0",           "azimut de rangee suppose nul"),
    (r"due[- ]south|plein sud|south-facing",
                                        "hypothese plein sud mentionnee en commentaire"),
    (r"min_beta_deg",                  "borne basse d'inclinaison propre au site"),
    (r"surface_azimuth|azimut de surface",
                                        "azimut de surface explicite (absent = implicite)"),
    (r"np\.sign\s*\(\s*lat",           "correction par signe de la latitude"),
    (r"Northern Hemisphere|hemisphere", "hemisphere mentionne"),
]

SCAN_DIRS = [REPO / "notebooks", REPO / "webapp" / "backend", REPO / "matlab"]


def map_assumptions():
    rows = []
    for d in SCAN_DIRS:
        if not d.exists():
            continue
        for f in sorted(d.rglob("*")):
            if f.suffix not in (".py", ".m") or not f.is_file():
                continue
            try:
                lines = f.read_text(encoding="utf-8", errors="replace").splitlines()
            except Exception:
                continue
            for i, line in enumerate(lines, 1):
                for pat, why in PATTERNS:
                    if re.search(pat, line, re.IGNORECASE):
                        rows.append({
                            "fichier": str(f.relative_to(REPO)),
                            "ligne": i,
                            "hypothese": why,
                            "extrait": line.strip()[:88],
                        })
    return pd.DataFrame(rows)


# ======================================================================
if __name__ == "__main__":
    pd.set_option("display.width", 190)
    pd.set_option("display.max_rows", 300)
    pd.set_option("display.max_colwidth", 90)

    print("=" * 96)
    print("C  Gel des references -- chemin de code actuel, 4 sites publies, scenario S0")
    print("=" * 96)
    frozen, cmp_df = freeze_references()
    Path("c_reference_frozen.json").write_text(json.dumps(frozen, indent=2, default=str))
    print(cmp_df.to_string(index=False))
    cmp_df.to_csv("c_reference_compare.csv", index=False)

    worst = cmp_df["ecart_pct"].abs().max()
    print()
    print(f"  Ecart max au CSV publie : {worst:.4f} %")
    print(f"  Reference gelee         : c_reference_frozen.json "
          f"({len(frozen)} sites, {sum(len(v['result']) + len(v['economics']) for v in frozen.values())} champs)")
    print(f"  Tolerance de non-regression retenue apres correction : {TOL_REL:g} (identite stricte,")
    print("  puisque les 4 sites sont dans l'hemisphere nord ou la correction est neutre).")

    print()
    print("=" * 96)
    print("D  Carte des hypotheses d'orientation dans le depot")
    print("=" * 96)
    amap = map_assumptions()
    if len(amap):
        print(amap.to_string(index=False))
    else:
        print("  (aucune occurrence)")
    amap.to_csv("d_assumption_map.csv", index=False)

    print()
    print("  Recapitulatif par hypothese :")
    for why, grp in amap.groupby("hypothese"):
        files = sorted(grp["fichier"].unique())
        print(f"    - {why}: {len(grp)} occurrence(s) dans {len(files)} fichier(s)")
