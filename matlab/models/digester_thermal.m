function th = digester_thermal(T_dig, T_amb, Q_heat_W, E_CH4_kWh_d, F, D, P)
%DIGESTER_THERMAL Lumped energy balance on the digester contents.
%
%   TH = DIGESTER_THERMAL(T_DIG, T_AMB, Q_HEAT_W, E_CH4_KWH_D, F, D, P)
%   returns
%
%     TH.dT_dt        digester temperature derivative (degC/d)
%     TH.Q_loss_W     conduction loss through the vessel wall (W)
%     TH.Q_flow_W     sensible heat to warm the influent (W)
%     TH.Q_reac_W     exothermic heat of digestion (W)
%     TH.C_th_kJ_K    thermal capacitance of the contents (kJ/K)
%     TH.tau_d        open-loop thermal time constant (d)
%
%   Energy balance:
%
%       C_th * dT/dt = Q_heat + Q_reac - Q_loss - Q_flow
%       C_th   = V * rho * Cp
%       Q_loss = U_eff * A_wall * (T_dig - T_amb)
%       Q_flow = Q_in * rho * Cp * (T_dig - T_in)
%       Q_reac = Q_reaction_frac * E_CH4
%
%   The wall area comes from the same cylindrical vessel (height = diameter)
%   the static model assumes, so the loss term is directly comparable.
%
%   What changes is that temperature is now a state rather than the proxy
%       T_dig_eff = clip(monthly mean T_amb + 5, 26, 37)
%   used in compute_AD(). That proxy has no memory: it says the digester at
%   Konya sits at 26.7 degC on annual average and cannot say whether it
%   holds 37 degC on a January night with the heater running. Here the
%   thermal capacitance -- of order 4 GJ/K for a 4 m3 vessel, a time
%   constant of days -- is exactly what buffers a heat supply that only
%   exists in daylight, and it is why the sizing question "how much PV
%   must be diverted to hold mesophilic conditions" has an answer at all.
%
%   See also DIGESTER_THERMOSTAT, AM2_GAS, APVAD_RHS.

W_TO_KJ_PER_DAY = 86.4;                     % 1 W sustained for 1 d = 86.4 kJ

C_th = F.V_dig_m3 * P.rho_substrate * P.Cp_substrate;      % kJ/K

U_eff  = P.U_eff_W_m2K * D.thermal.insulation_scale;
Q_loss = U_eff * F.A_wall_m2 * (T_dig - T_amb);            % W

% Sensible heating of the influent: kJ/d -> W
Q_flow_kJ_d = F.Q_in_m3_d * P.rho_substrate * P.Cp_substrate ...
              * max(T_dig - D.feed.T_in_degC, 0);
Q_flow = Q_flow_kJ_d / W_TO_KJ_PER_DAY;

% Exothermic heat of digestion, as a fraction of the methane energy
% produced -- the same closure the static model uses (Q_reaction_frac).
Q_reac = P.Q_reaction_frac * E_CH4_kWh_d * (1000/24);      % W

net_W = Q_heat_W + Q_reac - Q_loss - Q_flow;
dT_dt = net_W * W_TO_KJ_PER_DAY / C_th;                    % degC/d

% Open-loop time constant: C_th / (UA + rho*Cp*Q_in)
UA_eq_kJ_dK = U_eff*F.A_wall_m2*W_TO_KJ_PER_DAY ...
              + F.Q_in_m3_d*P.rho_substrate*P.Cp_substrate;

th.dT_dt     = dT_dt;
th.Q_loss_W  = Q_loss;
th.Q_flow_W  = Q_flow;
th.Q_reac_W  = Q_reac;
th.C_th_kJ_K = C_th;
th.tau_d     = C_th / max(UA_eq_kJ_dK, eps);
end
