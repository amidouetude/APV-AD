"""
Request isolation for arbitrary-point evaluation.

Guards the 2026-09 concurrency fix in model_bridge: a runtime site is
visible only inside the context that registered it, instead of being
written into the shared config_00.SITES dict.

TEST DESIGN -- read this before changing the parameters below.
The phase-0 audit measured that a concurrency test which varies the
COORDINATES detects nothing: api.py builds its key as
f"_runtime_{lat}_{lon}", so distinct cities give distinct keys and the
threads never contend. Under the old code that control case passed 8/8
while the real case failed 11/12. These tests therefore hammer the SAME
coordinates, which is also the realistic load pattern -- two users on the
same city, or one user double-clicking.
"""
import copy
import sys
import threading
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "notebooks"))
sys.path.insert(0, str(ROOT / "webapp" / "backend"))

import config_00  # noqa: E402
from model_bridge import (  # noqa: E402
    UnknownSiteError, evaluate_point, runtime_site,
)
import simulation_functions as sim  # noqa: E402

X = [16.37, 3.026, 1.501, 2.64, 31.8, 0.0502]
KEY = "_runtime_12.37_-1.53"          # exactly the key api.py builds
SAMPLE = ROOT / "data" / "ci_sample" / "hourly_Konya_sample.csv"


@pytest.fixture(scope="module")
def hourly():
    if not SAMPLE.exists():
        pytest.skip(f"CI sample not available at {SAMPLE}")
    df = pd.read_csv(SAMPLE)
    df["time_UTC"] = pd.to_datetime(df["time_UTC"])
    return df.set_index("time_UTC")


def _site(par_sat, yield_t_ha):
    s = copy.deepcopy(dict(config_00.SITES["Ouagadougou"]))
    s["PAR_sat"] = par_sat
    s["fruit_yield_t_ha"] = yield_t_ha
    return s


def _run_concurrently(jobs):
    """jobs = [(label, site_dict, par_sat)] on the SAME key."""
    results, errors = {}, {}
    barrier = threading.Barrier(len(jobs))

    def worker(label, site, par_sat):
        try:
            barrier.wait()
            results[label] = evaluate_point(KEY, site, worker.df, X,
                                             scenario="S0", PAR_sat=par_sat)
        except Exception as exc:                      # noqa: BLE001
            errors[label] = f"{type(exc).__name__}: {exc}"

    return results, errors, worker


# ======================================================================
def test_same_coordinates_different_crops_stay_isolated(hourly):
    """
    The case the old code failed: one key, two different site dicts. Every
    thread must get the result its OWN parameters imply.
    """
    configs = {"tomato": (174.0, 60.0), "maize": (400.0, 12.0)}

    serial = {name: evaluate_point(KEY, _site(*cfg), hourly, X,
                                    scenario="S0", PAR_sat=cfg[0])
              for name, cfg in configs.items()}
    assert serial["tomato"]["eLER_defA"] != serial["maize"]["eLER_defA"], (
        "the two configurations must differ, otherwise the test proves nothing"
    )

    jobs = [(f"{name}#{i}", _site(*cfg), cfg[0])
            for i in range(6) for name, cfg in configs.items()]
    results, errors, worker = _run_concurrently(jobs)
    worker.df = hourly
    threads = [threading.Thread(target=worker, args=j) for j in jobs]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert not errors, f"no thread may fail: {errors}"
    assert len(results) == len(jobs)

    for label, got in results.items():
        expected = serial[label.split("#")[0]]
        assert got["eLER_defA"] == pytest.approx(expected["eLER_defA"], abs=1e-12), (
            f"{label} got another request's parameters"
        )
        assert got["biogas_total_MWh_ha"] == pytest.approx(
            expected["biogas_total_MWh_ha"], abs=1e-12)


def test_same_coordinates_identical_params_all_succeed(hourly):
    """
    The other old failure: a thread's cleanup removed the key while its
    neighbours were still mid-pipeline, raising KeyError. Nobody may fail.
    """
    site = _site(174.0, 60.0)
    jobs = [(f"t{i}", copy.deepcopy(site), 174.0) for i in range(10)]
    results, errors, worker = _run_concurrently(jobs)
    worker.df = hourly
    threads = [threading.Thread(target=worker, args=j) for j in jobs]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert not errors, f"{len(errors)}/{len(jobs)} threads failed: {errors}"
    values = {round(r["eLER_defA"], 10) for r in results.values()}
    assert len(values) == 1, f"identical inputs gave different answers: {values}"


# ======================================================================
def test_runtime_site_does_not_touch_the_catalog(hourly):
    """The shared catalog must stay exactly the four published sites."""
    before = sorted(dict.keys(config_00.SITES))
    with runtime_site(KEY, _site(174.0, 60.0)):
        assert KEY in config_00.SITES            # visible inside
        assert KEY not in dict.keys(config_00.SITES)   # but not in the catalog
    assert KEY not in config_00.SITES            # gone outside
    assert sorted(dict.keys(config_00.SITES)) == before


def test_pipeline_and_config_share_one_registry():
    """
    simulation_functions.py does `from config_00 import SITES`, binding its
    own name. If only config_00.SITES were swapped, the pipeline would keep
    reading the old dict and the overlay would be silently ignored.
    """
    assert sim.SITES is config_00.SITES


def test_unknown_site_raises_a_readable_error():
    """Configuration failures must be explicit, not a bare KeyError."""
    with pytest.raises(UnknownSiteError) as exc:
        _ = config_00.SITES["_no_such_site"]
    assert "Unknown site" in str(exc.value)
    assert issubclass(UnknownSiteError, KeyError)   # backward compatible


def test_overlay_nesting_is_not_destructive(hourly):
    """An inner scope must not erase an outer one."""
    with runtime_site("outer", _site(174.0, 60.0)):
        with runtime_site("inner", _site(400.0, 12.0)):
            assert "outer" in config_00.SITES
            assert "inner" in config_00.SITES
        assert "outer" in config_00.SITES
        assert "inner" not in config_00.SITES


# ======================================================================
# Reserves soulevees en revue : le MECANISME, pas seulement le scenario
# mesure. Les 24 threads ci-dessus prouvent UN cas ; ces tests visent les
# proprietes generales dont dependront les taches de fond quand JOBS sera
# remplace par une file de travailleurs.
# ======================================================================
def test_background_thread_carries_its_own_context(hourly):
    """
    api.py lance ses jobs via threading.Thread. Un thread N'HERITE PAS du
    contexte de son createur : la surcouche doit donc etre posee DANS le
    thread, pas avant son demarrage. Ce test fige ce contrat -- si un jour
    quelqu'un enregistre le site dans le handler puis lance le thread, il
    doit echouer ici et non en production.
    """
    seen = {}

    def job(label, par_sat):
        # surcouche posee ici, dans le thread, comme le fait optimize_site
        with runtime_site(KEY, _site(par_sat, 60.0)):
            seen[label] = config_00.SITES[KEY]["PAR_sat"]

    threads = [threading.Thread(target=job, args=(f"j{i}", 100.0 + i))
               for i in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert seen == {f"j{i}": 100.0 + i for i in range(8)}
    assert KEY not in config_00.SITES          # rien n'a fui vers le thread principal


def test_outer_context_is_not_visible_to_a_child_thread():
    """
    Corollaire explicite : poser la surcouche AVANT de demarrer un thread ne
    la rend pas visible dedans. Documente pourquoi optimize_site et
    run_pareto ouvrent chacun leur propre scope.
    """
    visible = {}

    def peek():
        visible["seen"] = KEY in config_00.SITES

    with runtime_site(KEY, _site(174.0, 60.0)):
        assert KEY in config_00.SITES
        t = threading.Thread(target=peek)
        t.start()
        t.join()

    assert visible["seen"] is False, (
        "un thread enfant a herite du contexte -- le contrat a change, "
        "relire model_bridge et les appelants de runtime_site"
    )


def test_no_leak_between_sequential_requests(hourly):
    """
    Deux requetes successives sur la meme cle avec des parametres differents
    doivent donner chacune SON resultat : aucune reference mutable ne doit
    survivre a la premiere.
    """
    a = evaluate_point(KEY, _site(174.0, 60.0), hourly, X, scenario="S0", PAR_sat=174.0)
    b = evaluate_point(KEY, _site(400.0, 12.0), hourly, X, scenario="S0", PAR_sat=400.0)
    a2 = evaluate_point(KEY, _site(174.0, 60.0), hourly, X, scenario="S0", PAR_sat=174.0)

    assert a["eLER_defA"] != pytest.approx(b["eLER_defA"])
    assert a2["eLER_defA"] == pytest.approx(a["eLER_defA"], abs=1e-12), (
        "la troisieme requete a vu un residu de la deuxieme"
    )
    assert KEY not in config_00.SITES


def test_concurrent_asyncio_tasks_are_isolated(hourly):
    """
    Meme garantie cote asyncio : chaque tache porte son propre contexte, donc
    des endpoints async concurrents sur la meme cle ne se marchent pas dessus.
    """
    import asyncio

    async def one(par_sat):
        with runtime_site(KEY, _site(par_sat, 60.0)):
            await asyncio.sleep(0)             # force un point de bascule
            return config_00.SITES[KEY]["PAR_sat"]

    async def go():
        return await asyncio.gather(*(one(100.0 + i) for i in range(8)))

    assert asyncio.run(go()) == [100.0 + i for i in range(8)]
    assert KEY not in config_00.SITES


def test_every_lookup_goes_through_the_registry():
    """
    Le registre ne sert a rien si un module garde une reference au
    dictionnaire d'origine. Verifie que tout module resolvant un site pointe
    bien le meme objet.
    """
    import model_bridge
    registry = config_00.SITES
    assert isinstance(registry, model_bridge._SiteRegistry)
    assert sim.SITES is registry

    for name, mod in list(sys.modules.items()):
        if not name.startswith(("config_00", "simulation_functions", "model_bridge",
                                 "optimizer", "pareto", "economics", "api")):
            continue
        other = getattr(mod, "SITES", None)
        if other is not None:
            assert other is registry, f"{name}.SITES n'est pas le registre partage"
