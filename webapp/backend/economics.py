"""
4E economic/environmental analysis, extracted from
notebooks/04_metrics.ipynb (cells 2-5). The notebook itself is NOT
modified or imported here -- this is a standalone re-implementation of
the exact same formulas, generalized to take a `site` dict directly
(instead of a `SITES[name]` lookup) so it works for both the 4 catalog
sites and model_bridge.evaluate_point()'s runtime site dicts for
arbitrary locations.

If 04_metrics.ipynb's formulas ever change, this file needs the same
edit by hand -- there is no automated sync (unlike simulation_functions.py,
which 02_simulation.ipynb auto-exports). Kept deliberately close to the
notebook's own variable names so a side-by-side diff stays easy.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "notebooks"))
from config_00 import PARAMS, SCENARIOS  # noqa: E402  (unmodified)


def compute_capex(result: dict) -> dict:
    """Verbatim from 04_metrics.ipynb cell 2 (compute_capex)."""
    P = PARAMS
    capex_pv = result["n_modules"] * P["P_STC_Wp"] * P["capex_pv_usd_Wp"]
    capex_ad = result["V_dig_m3"] * P["capex_ad_usd_m3"]
    capex_agri = P["capex_agri_usd_ha"]
    capex_bop = P["bop_fraction"] * (capex_pv + capex_ad)
    capex_tot = capex_pv + capex_ad + capex_agri + capex_bop
    opex_yr = capex_tot * P["opex_fraction"]
    return dict(
        CAPEX_PV_kUSD=round(capex_pv / 1e3, 2),
        CAPEX_AD_kUSD=round(capex_ad / 1e3, 2),
        CAPEX_agri_kUSD=round(capex_agri / 1e3, 2),
        CAPEX_BOP_kUSD=round(capex_bop / 1e3, 2),
        CAPEX_total_kUSD=round(capex_tot / 1e3, 2),
        opex_kUSD_yr=round(opex_yr / 1e3, 2),
        _capex=capex_tot, _opex=opex_yr,
    )


def compute_revenue(site: dict, result: dict, scenario: str) -> dict:
    """Verbatim from 04_metrics.ipynb cell 3 (compute_revenue); `cfg=SITES[site]` becomes the `site` dict passed in."""
    P, cfg, sc = PARAMS, site, SCENARIOS[scenario]
    rev_pv = result["PV_sold_MWh_ha"] * 1000 * cfg["p_elec_usd_kWh"]
    fruit_t = cfg["fruit_yield_t_ha"] * result["LER_crop"]
    rev_tomato = fruit_t * cfg["p_tomato_usd_t"]
    E_bg_adj = result["E_biogas_kWh_ha"] * sc["BMP_mult"] * sc["eta_cap"]
    if sc["biogas_price"] == "biomethane":
        p_bg = cfg["p_elec_usd_kWh"] * P["biomethane_premium_factor"]
    elif sc["biogas_price"] == "LPG":
        p_bg = P["p_LPG_usd_kWh"]
    else:
        p_bg = cfg["p_biogas_usd_kWh"]
    rev_biogas = E_bg_adj * p_bg
    rev_digest = 0.0
    if sc["digestate_N"]:
        rev_digest = (fruit_t * 1000 * cfg["R_res"]
                      * P["N_content_kg_per_kg_residue"]
                      * P["N_price_usd_kg"])
    rev_carbon = 0.0
    if sc["carbon_credit"]:
        CO2 = (result["PV_sold_MWh_ha"] * cfg["grid_ef_tCO2_MWh"]
               + E_bg_adj / 1000 * P["ef_natural_gas_tCO2_MWh"])
        rev_carbon = CO2 * cfg["carbon_price_usd_tCO2"]
    gross = rev_pv + rev_tomato + rev_biogas + rev_digest + rev_carbon
    return dict(rev_pv=rev_pv, rev_tomato=rev_tomato, rev_biogas=rev_biogas,
                rev_digest=rev_digest, rev_carbon=rev_carbon, gross=gross,
                gross_kUSD=round(gross / 1e3, 3))


def _npv(capex, net_cf, r, N=None):
    """Verbatim from 04_metrics.ipynb cell 4."""
    N = N or PARAMS["project_life_yr"]
    if r == 0:
        return -capex + net_cf * N
    return -capex + net_cf * (1 - (1 + r) ** (-N)) / r


def _irr(capex, net_cf, N=None):
    """Verbatim from 04_metrics.ipynb cell 4 (bisection search)."""
    N = N or PARAMS["project_life_yr"]
    if net_cf <= 0 or _npv(capex, net_cf, 0, N) < 0:
        return float("nan")
    lo, hi = 1e-6, 10.0
    for _ in range(120):
        mid = (lo + hi) / 2
        if _npv(capex, net_cf, mid, N) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _lcoe(capex_pv, opex_pv, E_kWh, r, N=None):
    """Verbatim from 04_metrics.ipynb cell 4."""
    N = N or PARAMS["project_life_yr"]
    if E_kWh <= 0 or r <= 0:
        return float("nan")
    ann = (1 - (1 + r) ** (-N)) / r
    return (capex_pv + opex_pv * ann) / (E_kWh * ann)


def compute_4e(site: dict, result: dict, scenario: str = "S0") -> dict:
    """
    Full per-(site, scenario) 4E row -- generalized from 04_metrics.ipynb
    cell 5's per-row loop body to take a `site` dict directly (works for
    both SITES[name] and model_bridge's runtime site dicts) and a single
    `result` dict from eLER_objective(..., return_full=True).
    """
    P = PARAMS
    disc = site["discount_rate"]
    cap = compute_capex(result)

    CO2_PV = result["PV_sold_MWh_ha"] * site["grid_ef_tCO2_MWh"]
    CO2_bg = (result["E_biogas_kWh_ha"] / 1000) * P["ef_natural_gas_tCO2_MWh"]
    CO2_tot = CO2_PV + CO2_bg
    c_val = CO2_tot * site["carbon_price_usd_tCO2"]

    cpv = cap["CAPEX_PV_kUSD"] * 1e3
    pv_opex_share = cpv / cap["_capex"]
    lcoe_val = _lcoe(cpv, cap["_opex"] * pv_opex_share,
                      result["PV_sold_MWh_ha"] * 1000, disc)

    cad = cap["CAPEX_AD_kUSD"] * 1e3
    ad_share = cad / cap["_capex"]
    ann = (1 - (1 + disc) ** (-P["project_life_yr"])) / disc if disc > 0 else P["project_life_yr"]
    cbg = ((cad + cap["_opex"] * ad_share * ann) / (result["E_biogas_kWh_ha"] * ann)
           if result["E_biogas_kWh_ha"] > 0 else float("nan"))

    rv = compute_revenue(site, result, scenario)
    net = rv["gross"] - cap["_opex"]
    NPV = _npv(cap["_capex"], net, disc)
    IRR = _irr(cap["_capex"], net)
    PB = cap["_capex"] / net if net > 0 else float("inf")

    return {
        "CAPEX_PV_kUSD": cap["CAPEX_PV_kUSD"],
        "CAPEX_AD_kUSD": cap["CAPEX_AD_kUSD"],
        "CAPEX_agri_kUSD": cap["CAPEX_agri_kUSD"],
        "CAPEX_BOP_kUSD": cap["CAPEX_BOP_kUSD"],
        "CAPEX_total_kUSD": cap["CAPEX_total_kUSD"],
        "gross_revenue_kUSD_yr": round(rv["gross"] / 1e3, 2),
        "opex_kUSD_yr": round(cap["_opex"] / 1e3, 2),
        "net_CF_kUSD_yr": round(net / 1e3, 2),
        "NPV_kUSD": round(NPV / 1e3, 1),
        "IRR_pct": round(IRR * 100, 2) if IRR == IRR else None,
        "payback_yr": round(PB, 1) if PB < 1000 else None,
        "LCOE_PV_USD_kWh": round(lcoe_val, 4),
        "cost_biogas_USD_MWh": round(cbg, 2) if cbg == cbg else None,
        "CO2_avoided_PV_tCO2_ha": round(CO2_PV, 1),
        "CO2_avoided_biogas_tCO2_ha": round(CO2_bg, 2),
        "CO2_avoided_total_tCO2_ha": round(CO2_tot, 1),
        "carbon_value_USD_ha": round(c_val, 1),
    }
