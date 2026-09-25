"""
Bi-objective NSGA-II Pareto search (maximise eLER and LER_water jointly),
generalized from notebooks/07_pareto_analysis.py to run against an
arbitrary site dict + hourly DataFrame via the same runtime-injection
mechanism as model_bridge.py/optimizer.py, instead of requiring one of
the 4 catalog site names.

notebooks/07_pareto_analysis.py itself is not modified or imported --
this reimplements its EnergyWaterProblem class with the same objective
signs and the same DE_SETTINGS["bounds"] design space.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "notebooks"))
import config_00  # noqa: E402
import simulation_functions as sim  # noqa: E402

from model_bridge import runtime_site  # noqa: E402  (installs the site registry)
import numpy as np  # noqa: E402


def run_pareto(key: str, site_dict: dict, hourly_df, pop_size: int = 50,
                n_gen: int = 80, seed: int = 42, PAR_sat: float = 174.0,
                progress_cb=None):
    """
    Run the real NSGA-II search against an arbitrary site.

    progress_cb(generation, n_gen), if given, is called once per
    generation (pymoo callback) -- the hook for streaming job progress.

    Returns (eLER: np.ndarray, LER_water_pct: np.ndarray, X: np.ndarray)
    -- one row per non-dominated solution on the front.
    """
    from pymoo.core.problem import Problem
    from pymoo.core.callback import Callback
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.optimize import minimize

    bounds = config_00.DE_SETTINGS["bounds"]
    lb = np.array([b[0] for b in bounds])
    ub = np.array([b[1] for b in bounds])

    class EnergyWaterProblem(Problem):
        """pymoo minimizes by convention -- both objectives negated,
        exactly as in 07_pareto_analysis.py's EnergyWaterProblem."""

        def __init__(self):
            super().__init__(n_var=6, n_obj=2, n_constr=0, xl=lb, xu=ub)

        def _evaluate(self, X, out, *args, **kwargs):
            F = np.zeros((X.shape[0], 2))
            for i, x in enumerate(X):
                result = sim.eLER_objective(x, hourly_df, key, PAR_sat=PAR_sat, return_full=True)
                if result is None:
                    F[i, 0] = 10.0
                    F[i, 1] = 10.0
                else:
                    F[i, 0] = -result["eLER_defA"]
                    F[i, 1] = -result["LER_water_season_pct"]
            out["F"] = F

    class ProgressCallback(Callback):
        def notify(self, algorithm):
            if progress_cb:
                progress_cb(algorithm.n_gen, n_gen)

    # Per-request site scope -- see model_bridge's module docstring.
    with runtime_site(key, site_dict):
        problem = EnergyWaterProblem()
        algorithm = NSGA2(pop_size=pop_size)
        res = minimize(problem, algorithm, ("n_gen", n_gen), seed=seed,
                        callback=ProgressCallback(), verbose=False)

        eLER = -res.F[:, 0]
        LER_water = -res.F[:, 1]
        return eLER, LER_water, res.X
