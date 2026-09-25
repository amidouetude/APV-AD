function P = apvad_params()
%APVAD_PARAMS Site-independent physical constants shared with the static model.
%
%   P = APVAD_PARAMS() returns the constants ported from
%   notebooks/config_00.py (PARAMS). These are the values the static
%   six-model Python pipeline uses; keeping them identical here is what
%   makes the static-vs-dynamic comparison in COMPARE_WITH_STATIC valid.
%
%   Parameters that exist only in the dynamic model (AM2 kinetics, thermal
%   capacitances, crop and soil-water parameters) live in APVAD_DYN_PARAMS.
%
%   See also APVAD_SITES, APVAD_DYN_PARAMS, COMPARE_WITH_STATIC.

% ---- PV module: Jinko Tiger Neo 550 Wp bifacial -----------------------
P.P_STC_Wp        = 550.0;      % rated power at STC (W)
P.bifaciality     = 0.70;       % rear-to-front power ratio (phi)
P.beta_p          = -0.0035;    % power temperature coefficient (1/K)
P.eta_loss        = 0.12;       % system losses (wiring, soiling, mismatch)
P.rho_ground      = 0.25;       % ground albedo
P.phi_rear        = 0.12;       % rear irradiance as fraction of GHI-reflected

% ---- PV module / field geometry --------------------------------------
P.module_length_m = 2.278;      % m
P.module_width_m  = 1.134;      % m
P.row_length_m    = 100.0;      % m
P.land_area_m2    = 10000;      % 1 ha

% ---- Faiman thermal model (Faiman 2008) ------------------------------
P.U0              = 25.0;       % W/(m2 K) free-convection coefficient
P.U1              = 6.84;       % W s/(m3 K) wind-dependent coefficient

% ---- Crop / PAR -------------------------------------------------------
P.f_PAR           = 0.48;       % PAR fraction of GHI
P.alpha_shade     = 0.40;       % ET reduction coefficient under shade

% ---- Anaerobic digestion (static-model reference values) --------------
P.BMP_ref         = 300.0;      % NmL CH4 / g VS at 37 degC (Lallement 2023)
P.T_ref_degC      = 37.0;       % reference digester temperature
P.theta_arrhenius = 0.069;      % Arrhenius coefficient (1/degC) (Pilarski 2025)
P.T_dig_min_degC  = 26.0;       % minimum effective digester temperature
P.T_amb_boost_degC= 5.0;        % thermal boost above ambient (static proxy)
P.U_eff_W_m2K     = 0.8;        % digester wall heat transfer coefficient
P.CH4_fraction    = 0.60;       % static-model methane fraction (vol/vol)
P.LHV_CH4_MJ_m3   = 35.8;       % lower heating value of methane (MJ/Nm3)
P.eta_elec_heat   = 0.95;       % electrical resistance heating efficiency
P.rho_substrate   = 1020.0;     % kg/m3
P.Cp_substrate    = 4.2;        % kJ/(kg K)
P.T_inlet_degC    = 15.0;       % substrate inlet temperature (degC)
P.Q_reaction_frac = 0.05;       % exothermic heat as fraction of biogas energy

% ---- 4E / economics ---------------------------------------------------
P.ef_natural_gas_tCO2_MWh = 0.202;
P.project_life_yr   = 20;
P.capex_pv_usd_Wp   = 0.63;
P.capex_ad_usd_m3   = 800.0;
P.capex_agri_usd_ha = 2500.0;
P.bop_fraction      = 0.15;
P.opex_fraction     = 0.015;

% ---- Design-variable bounds (same order as the Python DE optimizer) ---
% [beta_deg, d_row_m, H_m_m, V_dig_m3, HRT_days, f_PV_heat]
P.design_lb = [15.0,  3.0, 1.5,  2.0, 20.0, 0.05];
P.design_ub = [40.0, 14.0, 5.0, 50.0, 40.0, 0.50];
P.design_names = {'beta_deg','d_row_m','H_m_m','V_dig_m3','HRT_days','f_PV_heat'};

% ---- Feasibility limits ----------------------------------------------
P.OLR_min = 1.5;                % kg VS/(m3 d)
P.OLR_max = 5.0;                % kg VS/(m3 d)
end
