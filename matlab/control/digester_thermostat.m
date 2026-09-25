function h = digester_thermostat(T_dig, P_array_W, E_CH4_kWh_d, F, D, P)
%DIGESTER_THERMOSTAT Proportional heating controller with PV-first dispatch.
%
%   H = DIGESTER_THERMOSTAT(T_DIG, P_ARRAY_W, E_CH4_KWH_D, F, D, P)
%   returns
%
%     H.Q_demand_W    heat called for by the controller (W)
%     H.Q_heat_W      heat actually delivered (W)
%     H.Q_pv_W        share supplied from the PV array (W)
%     H.Q_boiler_W    share supplied by the biogas boiler (W)
%     H.deficit_W     unmet demand (W)
%     H.P_pv_heat_W   PV power diverted to heating (W)
%     H.P_pv_sold_W   PV power exported (W)
%     H.biogas_used_kWh_d  methane energy burned for heat (kWh/d)
%
%   Control law -- proportional with saturation:
%       Q_demand = clip( Kp*(T_set - T_dig), 0, Q_max )
%
%   Dispatch order: photovoltaic electricity first, then the biogas boiler
%   as backup. The reason for that order is that PV heat is the whole point
%   of the coupling, and because f_PV_heat is a DESIGN variable in the
%   existing optimiser, this controller is where that variable finally does
%   something. In the static model f_PV_heat merely relabels a fraction of
%   annual PV energy from "sold" to "heat" -- an accounting split that
%   cannot fail. Here the diverted power has to arrive when the digester is
%   actually cold, and at night it does not arrive at all. The gap between
%   Q_demand and Q_heat is what the biogas boiler, or a heat store, has to
%   cover, and H.deficit_W records what neither covered.
%
%   The optimiser's preference for the lower bound f_PV_heat ~ 0.05 at all
%   four sites is exactly the kind of result this makes testable: it may be
%   correct, or it may be an artefact of an allocation that never has to
%   respect timing.
%
%   See also DIGESTER_THERMAL, SMOOTH_MIN, APVAD_RHS.

th = D.thermal;

% Heater capacity and gain are sized per design point in BUILD_FORCING.
Q_max = F.Q_max_W;
Kp    = F.Kp_W_per_K;

% Smoothing width for the saturations, as a fraction of heater capacity.
% This must scale with the quantity: a fixed 1 W width on a 1 kW heater is
% a kink in all but name, and the solver pays for it at every sunrise and
% sunset when PV availability crosses heat demand -- twice a day, every day
% of the year. Two percent of capacity is below any physical resolution the
% model claims and costs the integration nothing.
eps_Q = 0.02 * Q_max;

% ---- Demand -----------------------------------------------------------
err = th.T_set_degC - T_dig;
Q_demand = Kp * err;
Q_demand = Q_demand * smooth_step(err, th.T_db_degC);   % no cooling
Q_demand = smooth_min(Q_demand, Q_max, eps_Q);
Q_demand = max(Q_demand, 0);

% ---- Supply 1: photovoltaic electricity -------------------------------
P_pv_heat_avail = F.design.f_PV_heat * P_array_W;       % W diverted
Q_pv_avail = P_pv_heat_avail * P.eta_elec_heat;         % W of heat
Q_pv = smooth_min(Q_demand, Q_pv_avail, eps_Q);
Q_pv = max(Q_pv, 0);

% Any diverted PV power the digester does not need is exported instead of
% being dumped, so the allocation is a dispatch decision, not a loss.
P_pv_heat_used = Q_pv / max(P.eta_elec_heat, eps);
P_pv_sold = max(P_array_W - P_pv_heat_used, 0);

% ---- Supply 2: biogas boiler -----------------------------------------
Q_rem = max(Q_demand - Q_pv, 0);
if th.use_biogas_boiler
    % kWh/d -> W :  1 kWh/d = 1000/24 W
    Q_bio_avail = E_CH4_kWh_d * (1000/24) * th.boiler_eta;
    Q_boiler = smooth_min(Q_rem, Q_bio_avail, eps_Q);
    Q_boiler = max(Q_boiler, 0);
else
    Q_boiler = 0;
end
biogas_used_kWh_d = Q_boiler / max(th.boiler_eta, eps) * (24/1000);

% ---- Totals -----------------------------------------------------------
Q_heat = Q_pv + Q_boiler;

h.Q_demand_W  = Q_demand;
h.Q_heat_W    = Q_heat;
h.Q_pv_W      = Q_pv;
h.Q_boiler_W  = Q_boiler;
h.deficit_W   = max(Q_demand - Q_heat, 0);
h.P_pv_heat_W = P_pv_heat_used;
h.P_pv_sold_W = P_pv_sold;
h.biogas_used_kWh_d = biogas_used_kWh_d;
end
