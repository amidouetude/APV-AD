"""
Real DE + L-BFGS-B optimizer, generalized from notebooks/03_optimizer.ipynb
(cells 2-3: same scipy.optimize.differential_evolution call, same
hyperparameters from config_00.DE_SETTINGS, same manual L-BFGS-B polish
step) to run against an arbitrary site dict + hourly DataFrame via
model_bridge's runtime-injection mechanism, instead of requiring one of
the 4 catalog site names baked into config_00.SITES.

notebooks/03_optimizer.ipynb itself is not modified or imported.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "notebooks"))
import config_00  # noqa: E402
import simulation_functions as sim  # noqa: E402
from scipy.optimize import differential_evolution, minimize  # noqa: E402


def optimize_site(key: str, site_dict: dict, hourly_df, scenario: str = "S0",
                   PAR_sat: float = 174.0, progress_cb=None, max_generations=None):
    """
    Run the real DE + L-BFGS-B optimizer against an arbitrary site.

    progress_cb(generation, best_eLER), if given, is called once per DE
    generation -- the hook for streaming a live convergence chart.

    Returns (x_opt: list[float], full_result: dict, history: list[float]).
    """
    config_00.SITES[key] = site_dict
    n_gen = max_generations or config_00.DE_SETTINGS["max_generations"]
    history = []

    def neg_eLER(x):
        val = sim.eLER_objective(x, hourly_df, key, scenario=scenario, PAR_sat=PAR_sat)
        return -val

    def callback(xk, convergence):
        val = -neg_eLER(xk)
        history.append(val)
        if progress_cb:
            progress_cb(len(history), val)
        return False

    try:
        DE = config_00.DE_SETTINGS
        result = differential_evolution(
            neg_eLER,
            bounds=DE["bounds"],
            maxiter=n_gen,
            popsize=DE["population_size"] // len(DE["bounds"]),
            mutation=DE["mutation_F"],
            recombination=DE["crossover_Cr"],
            seed=DE["seed"],
            callback=callback,
            tol=1e-6,
            polish=False,   # manual L-BFGS-B polish below, matching the notebook
            workers=1,
        )

        polish_result = minimize(
            neg_eLER, x0=result.x, method="L-BFGS-B",
            bounds=[(lo, hi) for lo, hi in DE["bounds"]],
            options={"maxiter": 500, "ftol": 1e-10},
        )
        x_opt = polish_result.x if polish_result.fun < result.fun else result.x

        full = sim.eLER_objective(x_opt, hourly_df, key, scenario=scenario,
                                   PAR_sat=PAR_sat, return_full=True)
        return list(x_opt), full, history
    finally:
        config_00.SITES.pop(key, None)
