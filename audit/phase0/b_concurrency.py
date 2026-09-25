"""
AUDIT PHASE 0 -- B : isolation des requetes concurrentes
=======================================================

LECTURE SEULE. Aucun fichier du depot n'est ecrit.

Ce que dit le code
------------------
model_bridge.evaluate_point() fait :

    register_runtime_site(key, site_dict)   # config_00.SITES[key] = site_dict
    try:
        return sim.eLER_objective(x, df, site=key, ...)
    finally:
        config_00.SITES.pop(key, None)

et api.py construit la cle ainsi :  key = f"_runtime_{lat}_{lon}"

La cle ne depend donc QUE des coordonnees. Or le site_dict depend aussi
de la culture (PAR_sat) et du site d'emprunt (borrow_from : rendement,
prix, taux d'actualisation, facteur d'emission...). Deux requetes sur le
MEME point avec une culture differente partagent la cle mais portent des
parametres differents.

Trois scenarios testes
----------------------
  B1  memes coordonnees, cultures differentes
      -> contamination silencieuse attendue (mauvais resultat, sans erreur)
  B2  coordonnees et parametres identiques
      -> le pop() d'un thread peut retirer la cle sous un autre : KeyError
  B3  coordonnees differentes
      -> devrait etre propre ; on verifie que c'est bien le cas
"""
import copy
import sys
import threading
import traceback
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent / "repo"
sys.path.insert(0, str(REPO / "notebooks"))
sys.path.insert(0, str(REPO / "webapp" / "backend"))

import config_00                      # noqa: E402
import simulation_functions as sim    # noqa: E402
from model_bridge import evaluate_point  # noqa: E402

X = [16.37, 3.026, 1.501, 2.64, 31.8, 0.0502]   # design d'Ouagadougou


def load_hourly(site):
    df = pd.read_csv(REPO / "outputs" / "csv" / f"hourly_{site}.csv")
    df["time_UTC"] = pd.to_datetime(df["time_UTC"])
    return df.set_index("time_UTC")


def make_site(par_sat, yield_t_ha):
    s = copy.deepcopy(config_00.SITES["Ouagadougou"])
    s["PAR_sat"] = par_sat
    s["fruit_yield_t_ha"] = yield_t_ha
    return s


# ----------------------------------------------------------------------
def run_threads(jobs):
    """jobs = [(label, key, site_dict, df, par_sat)] -> {label: result|exception}"""
    out, errs = {}, {}
    barrier = threading.Barrier(len(jobs))

    def worker(label, key, site, df, par_sat):
        try:
            barrier.wait()                       # depart simultane
            out[label] = evaluate_point(key, site, df, X, scenario="S0",
                                         PAR_sat=par_sat)
        except Exception as exc:
            errs[label] = f"{type(exc).__name__}: {exc}"

    ts = [threading.Thread(target=worker, args=j) for j in jobs]
    for t in ts:
        t.start()
    for t in ts:
        t.join()
    return out, errs


# ======================================================================
def b1_same_coords_different_crop(df, n_pairs=6):
    """Contamination silencieuse : meme cle, site_dict differents."""
    key = "_runtime_12.37_-1.53"          # exactement la cle que forge api.py

    cfg_a = ("tomate",  174.0, 60.0)
    cfg_b = ("mais",    400.0, 12.0)

    # references en serie, chacune seule
    ref = {}
    for label, par, yld in (cfg_a, cfg_b):
        ref[label] = evaluate_point(key, make_site(par, yld), df, X,
                                     scenario="S0", PAR_sat=par)

    jobs = []
    for i in range(n_pairs):
        for label, par, yld in (cfg_a, cfg_b):
            jobs.append((f"{label}#{i}", key, make_site(par, yld), df, par))

    out, errs = run_threads(jobs)

    contaminated, clean = [], []
    for label, res in out.items():
        base = label.split("#")[0]
        exp = ref[base]
        same = (res is not None and exp is not None
                and abs(res["eLER_defA"] - exp["eLER_defA"]) < 1e-9
                and abs(res["biogas_total_MWh_ha"] - exp["biogas_total_MWh_ha"]) < 1e-9)
        (clean if same else contaminated).append(label)

    return {
        "ref_eLER": {k: round(v["eLER_defA"], 4) for k, v in ref.items()},
        "ref_biogas": {k: round(v["biogas_total_MWh_ha"], 3) for k, v in ref.items()},
        "n_threads": len(jobs), "n_ok": len(clean),
        "n_contamines": len(contaminated), "contamines": sorted(contaminated)[:8],
        "erreurs": errs,
    }


def b2_same_key_identical(df, n=10):
    """pop() premature : un thread retire la cle pendant qu'un autre calcule."""
    key = "_runtime_12.37_-1.53"
    site = make_site(174.0, 60.0)
    jobs = [(f"t{i}", key, copy.deepcopy(site), df, 174.0) for i in range(n)]
    out, errs = run_threads(jobs)
    return {"n_threads": n, "n_reussis": len(out), "n_erreurs": len(errs),
            "erreurs": dict(list(errs.items())[:5])}


def b3_different_coords(df, n=8):
    """Coordonnees distinctes -> cles distinctes : devrait etre propre."""
    jobs, ref = [], {}
    for i in range(n):
        lat = 12.37 + i * 0.5
        key = f"_runtime_{lat}_-1.53"
        site = make_site(174.0, 60.0)
        site["lat"] = lat
        ref[f"p{i}"] = evaluate_point(key, copy.deepcopy(site), df, X,
                                       scenario="S0", PAR_sat=174.0)
        jobs.append((f"p{i}", key, copy.deepcopy(site), df, 174.0))

    out, errs = run_threads(jobs)
    bad = [k for k, v in out.items()
           if abs(v["eLER_defA"] - ref[k]["eLER_defA"]) > 1e-9]
    return {"n_threads": n, "n_ok": len(out) - len(bad),
            "n_incoherents": len(bad), "erreurs": errs}


# ======================================================================
if __name__ == "__main__":
    df = load_hourly("Ouagadougou")
    print(f"Jeu de donnees : Ouagadougou, {len(df)} heures\n")

    print("=" * 74)
    print("B1  Memes coordonnees, cultures differentes")
    print("=" * 74)
    r1 = b1_same_coords_different_crop(df)
    print(f"  references en serie  : eLER {r1['ref_eLER']}")
    print(f"                         biogaz {r1['ref_biogas']} MWh/ha")
    print(f"  threads              : {r1['n_threads']}")
    print(f"  resultats corrects   : {r1['n_ok']}")
    print(f"  resultats CONTAMINES : {r1['n_contamines']}  {r1['contamines']}")
    if r1["erreurs"]:
        print(f"  erreurs              : {r1['erreurs']}")
    verdict1 = "CONFIRME" if r1["n_contamines"] else "non reproduit"
    print(f"  -> contamination silencieuse : {verdict1}")

    print()
    print("=" * 74)
    print("B2  Memes coordonnees, parametres identiques (pop premature)")
    print("=" * 74)
    r2 = b2_same_key_identical(df)
    print(f"  threads {r2['n_threads']} | reussis {r2['n_reussis']} | erreurs {r2['n_erreurs']}")
    for k, v in r2["erreurs"].items():
        print(f"    {k}: {v}")
    verdict2 = "CONFIRME" if r2["n_erreurs"] else "non reproduit"
    print(f"  -> KeyError sous concurrence : {verdict2}")

    print()
    print("=" * 74)
    print("B3  Coordonnees differentes (temoin)")
    print("=" * 74)
    r3 = b3_different_coords(df)
    print(f"  threads {r3['n_threads']} | coherents {r3['n_ok']} | incoherents {r3['n_incoherents']}")
    if r3["erreurs"]:
        print(f"  erreurs : {r3['erreurs']}")
    print(f"  -> isolation par coordonnees distinctes : "
          f"{'OK' if not r3['n_incoherents'] and not r3['erreurs'] else 'DEFAILLANTE'}")

    import json
    Path("b_concurrency_results.json").write_text(
        json.dumps({"B1": r1, "B2": r2, "B3": r3}, indent=2, default=str))
