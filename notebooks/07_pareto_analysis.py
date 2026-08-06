"""
07_pareto_analysis.py
=============================================================================
Generates Figure 10 (eLER vs. LER_water Pareto front) referenced in
paper/sections/05_discussion.tex ("Beyond the scalar optimum: a Pareto
perspective").

Implements the bi-objective NSGA-II search already specified in the
manuscript's Methodology (Section 3.9.3, "Multi-Objective Extension") but
not previously run: maximises eLER and LER_water jointly, using the exact
same six-dimensional design space and the same simulation_functions.py
used for the primary scalar Differential Evolution optimization in
03_optimizer.ipynb.

Requires: pymoo (pip install pymoo), matplotlib, and the standard
project dependencies (numpy, pandas). Run AFTER 01_load_data.ipynb has
produced outputs/csv/hourly_{site}.csv, and after 02_simulation.ipynb has
been run at least once (so notebooks/simulation_functions.py exists).

Usage:
    python 07_pareto_analysis.py [--site Konya] [--pop 50] [--gen 80]

Output:
    outputs/figures/fig10_pareto_eler_water.png
=============================================================================
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── Import project modules exactly as 03_optimizer.ipynb does ──────────────
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from config_00 import SITES, DE_SETTINGS, CSV_DIR, FIG_DIR
from simulation_functions import eLER_objective


def load_site_data(site: str) -> pd.DataFrame:
    """Loads the pre-processed hourly climate CSV the same way
    02_simulation.ipynb and 03_optimizer.ipynb do — NOT raw PVGIS data.
    Requires 01_load_data.ipynb to have been run first."""
    path = CSV_DIR / f"hourly_{site}.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run 01_load_data.ipynb first to produce "
            f"the pre-processed hourly climate CSVs."
        )
    df = pd.read_csv(path, index_col=0, parse_dates=True)
    return df


def run_pareto_search(site: str, pop_size: int, n_gen: int, seed: int = 42):
    """Runs the bi-objective NSGA-II search (maximise eLER, maximise
    LER_water) over the same six-dimensional design space and bounds used
    for the primary scalar DE optimization (DE_SETTINGS['bounds'])."""
    from pymoo.core.problem import Problem
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.optimize import minimize

    df = load_site_data(site)
    bounds = DE_SETTINGS["bounds"]
    lb = np.array([b[0] for b in bounds])
    ub = np.array([b[1] for b in bounds])

    class EnergyWaterProblem(Problem):
        """pymoo minimises by convention, so both objectives are negated
        (eLER_objective and LER_water_season_pct are both maximised in the
        original sense)."""

        def __init__(self):
            super().__init__(n_var=6, n_obj=2, n_constr=0, xl=lb, xu=ub)

        def _evaluate(self, X, out, *args, **kwargs):
            F = np.zeros((X.shape[0], 2))
            for i, x in enumerate(X):
                result = eLER_objective(x, df, site, return_full=True)
                if result is None:
                    # Infeasible design (e.g. OLR or tilt constraint
                    # violated) — penalize both objectives.
                    F[i, 0] = 10.0
                    F[i, 1] = 10.0
                else:
                    F[i, 0] = -result["eLER_defA"]
                    F[i, 1] = -result["LER_water_season_pct"]
            out["F"] = F

    problem = EnergyWaterProblem()
    algorithm = NSGA2(pop_size=pop_size)
    res = minimize(problem, algorithm, ("n_gen", n_gen), seed=seed, verbose=True)

    eLER = -res.F[:, 0]
    LER_water = -res.F[:, 1]
    return eLER, LER_water, res.X


def plot_pareto_front(eLER, LER_water, site, scalar_eLER_optimum, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = np.argsort(eLER)
    eLER_s, LER_water_s = eLER[order], LER_water[order]

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.scatter(eLER_s, LER_water_s, c="#2a6f97", s=28, zorder=3,
               edgecolor="white", linewidth=0.5)
    ax.plot(eLER_s, LER_water_s, "--", color="#2a6f97", alpha=0.4,
            zorder=2, linewidth=1)

    ax.axvline(scalar_eLER_optimum, color="grey", linestyle=":",
               linewidth=1, zorder=1)
    y_annot = np.percentile(LER_water_s, 40)
    ax.annotate(
        f"Scalar eLER optimum\n(Table 9, {site})",
        xy=(scalar_eLER_optimum, y_annot),
        xytext=(eLER_s.min() + 0.10 * (eLER_s.max() - eLER_s.min()), y_annot + 4),
        fontsize=8.5, color="dimgrey",
        arrowprops=dict(arrowstyle="->", color="dimgrey", lw=0.8),
    )

    ax.set_xlabel("eLER (Definition A)", fontsize=11)
    ax.set_ylabel(r"$LER_{\mathrm{water}}$ season (%)", fontsize=11)
    ax.set_title(
        f"Pareto Front: eLER vs. Water-Saving Trade-off\n"
        f"({site}, NSGA-II, {len(eLER)} solutions)",
        fontsize=10,
    )
    ax.grid(alpha=0.25, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    print(f"Saved figure -> {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site", default="Konya", choices=list(SITES.keys()))
    parser.add_argument("--pop", type=int, default=50, help="NSGA-II population size")
    parser.add_argument("--gen", type=int, default=80, help="NSGA-II generations")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--scalar-eler-optimum", type=float, default=None,
        help="eLER value of the scalar DE optimum to annotate (Table 9). "
             "If omitted, uses this run's own maximum eLER on the front.",
    )
    parser.add_argument(
        "--out", type=Path, default=None,
        help="Output PNG path (default: outputs/figures/fig10_pareto_eler_water.png)",
    )
    args = parser.parse_args()

    print(f"Running NSGA-II for {args.site}: pop={args.pop}, gen={args.gen} "
          f"({args.pop * args.gen} total evaluations)...")
    eLER, LER_water, X = run_pareto_search(args.site, args.pop, args.gen, args.seed)

    scalar_opt = args.scalar_eler_optimum
    if scalar_opt is None:
        scalar_opt = eLER.max()
        print(f"No --scalar-eler-optimum given; using this run's own "
              f"max eLER ({scalar_opt:.3f}) for the reference line. "
              f"Pass the Table 9 value for {args.site} explicitly to "
              f"match the manuscript exactly.")

    out_path = args.out or (FIG_DIR / "fig10_pareto_eler_water.png")
    plot_pareto_front(eLER, LER_water, args.site, scalar_opt, out_path)

    print(f"\neLER range on front:      {eLER.min():.3f} -- {eLER.max():.3f}")
    print(f"LER_water range on front: {LER_water.min():.2f}% -- {LER_water.max():.2f}%")


if __name__ == "__main__":
    main()
