function [dx, out] = apvad_rhs(t, x, F, D, P)
%APVAD_RHS Right-hand side of the coupled APV-AD dynamic model.
%
%   DX = APVAD_RHS(T, X, F, D, P) returns the state derivative for ODE15S.
%   [DX, OUT] = APVAD_RHS(...) additionally returns a struct of algebraic
%   quantities (powers, flows, pH, stress factors) for postprocessing. The
%   second output is never requested inside the integration loop.
%
%   Inputs
%     T   time in days since 1 January 00:00
%     X   state vector, see APVAD_STATES
%     F   forcing struct from BUILD_FORCING
%     D   dynamic parameters, APVAD_DYN_PARAMS
%     P   shared physical constants, APVAD_PARAMS
%
%   Structure of the coupling
%   -------------------------
%   Climate -> solar geometry and shading (precomputed in F)
%       |                    |
%       |                    +--> PAR beneath the array -> crop growth
%       |                    +--> ET0 * (1 - alpha*F_shad) -> soil water
%       |                                                       |
%       +--> POA -> PV module (thermal state) -> array power     v
%                                    |                        Ks -> growth
%                                    +--> thermostat --> digester heating
%                                    +--> export                    |
%                                                                   v
%   crop harvest -> feedstock store -> AM2 digester <-- temperature state
%                                          |
%                                          +--> methane -> boiler backup
%
%   The two feedback loops that do not exist in the static model are worth
%   naming, because they are the reason for doing this at all:
%
%   1. Thermal-biological loop. Digester temperature sets the growth rates
%      through GROWTH_TEMP_FACTOR; the resulting methane both heats the
%      vessel (exothermic term and boiler backup) and is the output being
%      measured. Cold digester -> slower methanogenesis -> less gas -> less
%      backup heat available. The static model breaks this loop by fixing
%      temperature from an ambient proxy.
%
%   2. Shade-water-growth loop. Shading cuts evapotranspiration, which
%      lowers root-zone depletion, which raises Ks, which raises growth,
%      which raises residue yield, which raises digester feed. The static
%      model computes the water saving but never lets it act on yield --
%      LER_crop there is a pure light ratio.
%
%   See also APVAD_STATES, RUN_DYNAMIC_SIM, BUILD_FORCING.

ix = apvad_states();

% =====================================================================
% 1. Unpack states (floored where a negative value is unphysical)
% =====================================================================
S1     = max(x(ix.S1), 0);
X1     = max(x(ix.X1), 0);
S2     = max(x(ix.S2), 0);
X2     = max(x(ix.X2), 0);
Z      = max(x(ix.Z),  0);
Cinorg = max(x(ix.C),  0);
T_dig  = x(ix.T_dig);
T_cell = x(ix.T_cell);
Dr     = x(ix.Dr);
GDD    = max(x(ix.GDD), 0);
B_dm   = max(x(ix.B_dm), 0);
M_res  = max(x(ix.M_res), 0);

% =====================================================================
% 2. Exogenous forcing at time t
% =====================================================================
% Time is wrapped onto the climate record so that a multi-year run repeats
% the same typical meteorological year and converges to a periodic steady
% state. The first year is spin-up and is discarded by RUN_DYNAMIC_SIM.
tq = mod(t, F.t_period);

T_amb    = F.T_amb(tq);
WS       = F.WS(tq);
G_eff    = F.G_eff(tq);
F_shad   = F.F_shad(tq);
PAR_AV   = F.PAR_AV(tq);
ET0      = F.ET0(tq);
rain     = F.rain(tq);
Kc       = F.Kc(tq);
inSeason = min(max(F.in_season(tq), 0), 1);
T_day    = F.T_day(tq);
hPulse   = F.harvest(tq);

% =====================================================================
% 3. Photovoltaic array
% =====================================================================
if ~D.pv.dynamic
    T_cell = T_amb + G_eff / (P.U0 + P.U1*max(WS, 0.5));   % algebraic Faiman
end
pv = pv_module(G_eff, T_amb, WS, T_cell, F.n_modules, D, P);

% =====================================================================
% 4. Digester biology and gas phase
% =====================================================================
kin = am2_kinetics(S1, S2, T_dig, D, P);
gas = am2_gas(S2, Z, Cinorg, kin.mu2, X2, F.V_dig_m3, D, P);

% =====================================================================
% 5. Heating controller and digester thermal balance
% =====================================================================
heat = digester_thermostat(T_dig, pv.P_array_W, gas.E_CH4_kWh_d, F, D, P);
therm = digester_thermal(T_dig, T_amb, heat.Q_heat_W, gas.E_CH4_kWh_d, ...
                         F, D, P);

% =====================================================================
% 6. Soil water and crop
% =====================================================================
wat  = soil_water(Dr, ET0, Kc, F_shad, rain, inSeason, D, P);
crop = crop_growth(GDD, B_dm, PAR_AV, F.site.PAR_sat, T_day, inSeason, ...
                   wat.Ks, hPulse, D);

% Harvested dry matter splits into marketable fruit and residues; only the
% residues -- at the site's VS/DM ratio -- reach the digester store.
residue_DM_kg_d = crop.harvest_kgDM_d * (1 - D.crop.HI);
residue_VS_kg_d = residue_DM_kg_d * F.site.VS_fraction;

% =====================================================================
% 7. Feedstock store and influent composition
% =====================================================================
fd = feed_scheduler(M_res, residue_VS_kg_d, F, D);

% =====================================================================
% 8. AM2 mass balances
% =====================================================================
a   = D.am2;
Dil = F.Q_in_m3_d / max(F.V_dig_m3, eps);        % dilution rate, 1/d

dS1 = Dil*(fd.S1in - S1) - a.k1*kin.mu1*X1;
dX1 = (kin.mu1 - a.alpha*Dil - a.kd1)*X1;
dS2 = Dil*(fd.S2in - S2) + a.k2*kin.mu1*X1 - a.k3*kin.mu2*X2;
dX2 = (kin.mu2 - a.alpha*Dil - a.kd2)*X2;
dZ  = Dil*(fd.Zin - Z);
dC  = Dil*(fd.Cin - Cinorg) - gas.qC_mmol_L_d ...
      + a.k4*kin.mu1*X1 + a.k5*kin.mu2*X2;

% =====================================================================
% 9. Assemble the derivative
% =====================================================================
dx = zeros(ix.n, 1);

dx(ix.S1)    = dS1;
dx(ix.X1)    = dX1;
dx(ix.S2)    = dS2;
dx(ix.X2)    = dX2;
dx(ix.Z)     = dZ;
dx(ix.C)     = dC;

dx(ix.T_dig) = therm.dT_dt;
dx(ix.T_cell)= pv.dTcell_dt;

dx(ix.Dr)    = wat.dDr_dt;
dx(ix.GDD)   = crop.dGDD_dt;
dx(ix.B_dm)  = crop.dB_dt;
dx(ix.M_res) = fd.dM_res_dt;

% Accumulators. 1 W sustained for one day is 24/1000 kWh.
dx(ix.E_pv)       = pv.P_array_W   * 24/1000;
dx(ix.E_ch4)      = gas.E_CH4_kWh_d;
dx(ix.E_heat)     = heat.Q_heat_W  * 24/1000;
dx(ix.E_pv_sold)  = heat.P_pv_sold_W * 24/1000;
dx(ix.E_ch4_heat) = heat.biogas_used_kWh_d;
dx(ix.W_irr)      = wat.irr_mm_d;
dx(ix.M_harv)     = crop.harvest_kgDM_d;

% =====================================================================
% 10. Algebraic outputs (postprocessing only)
% =====================================================================
if nargout > 1
    out.t = t;

    out.T_amb = T_amb;  out.WS = WS;  out.G_eff = G_eff;
    out.F_shad = F_shad; out.PAR_AV = PAR_AV; out.ET0 = ET0;
    out.Kc = Kc; out.in_season = inSeason; out.T_day = T_day;

    out.T_cell     = T_cell;
    out.T_cell_ss  = pv.T_cell_ss;
    out.P_array_W  = pv.P_array_W;
    out.P_pv_sold_W= heat.P_pv_sold_W;
    out.P_pv_heat_W= heat.P_pv_heat_W;

    out.mu1 = kin.mu1;  out.mu2 = kin.mu2;  out.fT = kin.fT;
    out.Q_CH4_Nm3_d = gas.Q_CH4_Nm3_d;
    out.Q_biogas_Nm3_d = gas.Q_biogas_Nm3_d;
    out.x_CH4 = gas.x_CH4;
    out.pH = gas.pH;
    out.E_CH4_kWh_d = gas.E_CH4_kWh_d;

    out.Q_demand_W = heat.Q_demand_W;
    out.Q_heat_W   = heat.Q_heat_W;
    out.Q_pv_W     = heat.Q_pv_W;
    out.Q_boiler_W = heat.Q_boiler_W;
    out.deficit_W  = heat.deficit_W;
    out.Q_loss_W   = therm.Q_loss_W;
    out.Q_flow_W   = therm.Q_flow_W;
    out.tau_dig_d  = therm.tau_d;

    out.Ks        = wat.Ks;
    out.ETc_mm_d  = wat.ETc_mm_d;
    out.ETc_open  = wat.ETc_open;
    out.irr_mm_d  = wat.irr_mm_d;
    out.TAW       = wat.TAW;
    out.RAW       = wat.RAW;

    out.LAI       = crop.LAI;
    out.f_i       = crop.f_i;
    out.APAR_MJ   = crop.APAR_MJ;
    out.growth_kgDM_d = crop.growth_kgDM_d;

    out.OLR_kgVS_m3d = fd.OLR_kgVS_m3d;
    out.VS_feed_kg_d = fd.VS_feed_kg_d;
    out.S1in         = fd.S1in;
    out.store_avail  = fd.avail;
    out.Dil          = Dil;
end
end
