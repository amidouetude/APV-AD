function pv = pv_module(G_eff, T_amb, WS, T_cell, n_modules, D, P)
%PV_MODULE Bifacial PV array power and module thermal dynamics.
%
%   PV = PV_MODULE(G_EFF, T_AMB, WS, T_CELL, N_MODULES, D, P) returns
%
%     PV.T_cell_ss    steady-state (Faiman) module temperature (degC)
%     PV.dTcell_dt    module temperature derivative (degC/d)
%     PV.P_array_W    array electrical power (W)
%     PV.P_module_W   power of a single module (W)
%     PV.eta_temp     temperature derating factor (-)
%
%   Steady-state module temperature (Faiman 2008):
%       T_ss = T_amb + G_eff / (U0 + U1*WS)
%
%   Dynamic module temperature -- a lumped first-order energy balance on
%   the module areal heat capacity:
%       C_areal * dT/dt = G_eff - (U0 + U1*WS)*(T_cell - T_amb) [- P_elec/A]
%
%   By construction the steady state of this ODE is exactly the Faiman
%   expression, so setting D.pv.dynamic = false (algebraic Faiman) recovers
%   the static Python model bit for bit. That equivalence is what
%   TESTPVTHERMALREDUCESTOFAIMAN checks, and it is the reason the PV
%   sub-model can be trusted as a control rather than a confound when the
%   dynamic and static annual yields are compared.
%
%   With C_areal ~ 18 kJ/(m2 K) and U ~ 25 W/(m2 K) the module time
%   constant is about 12 minutes. Against hourly forcing that lag changes
%   annual yield by well under a percent; it matters for sub-hourly studies
%   and for anything that looks at ramp rates, which is precisely what a
%   dynamic model is for.
%
%   Array power:
%       P = P_STC * N * (1 - eta_loss) * [1 + beta_p*(T_cell - 25)]
%           * G_eff/1000
%
%   Reference
%     Faiman, D. (2008). Prog. Photovolt: Res. Appl., 16(4), 307-315.
%
%   See also APVAD_RHS, POA_IRRADIANCE.

U = P.U0 + P.U1 * max(WS, 0.5);          % W/(m2 K), wind floor as in Python

T_ss = T_amb + G_eff ./ U;

% ---- Electrical power at the current module temperature --------------
eta_temp   = max(1 + P.beta_p*(T_cell - 25), 0);
P_module_W = max(P.P_STC_Wp * (1 - P.eta_loss) .* eta_temp .* G_eff/1000, 0);
P_array_W  = P_module_W * n_modules;

% ---- Thermal derivative (degC/d) -------------------------------------
% 1 W/m2 sustained for one day delivers 86.4 kJ/m2.
if D.pv.dynamic
    q_net = G_eff - U.*(T_cell - T_amb);                    % W/m2
    if D.pv.include_Pelec
        A_mod = P.module_length_m * P.module_width_m;
        q_net = q_net - P_module_W / A_mod;
    end
    dTcell_dt = 86.4 * q_net / D.pv.C_areal_kJ_m2K;         % degC/d
else
    dTcell_dt = 0;                                          % held algebraic
end

pv.T_cell_ss  = T_ss;
pv.dTcell_dt  = dTcell_dt;
pv.P_array_W  = P_array_W;
pv.P_module_W = P_module_W;
pv.eta_temp   = eta_temp;
pv.U_W_m2K    = U;
end
