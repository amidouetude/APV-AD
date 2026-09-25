function D = apvad_dyn_params()
%APVAD_DYN_PARAMS Parameters that exist only in the dynamic APV-AD model.
%
%   D = APVAD_DYN_PARAMS() returns kinetic, thermal, crop and soil-water
%   parameters required by the continuous-time model. None of these appear
%   in the static Python pipeline, which represents the digester as an
%   annual lumped energy balance and the crop as a season-integrated PAR
%   ratio.
%
%   IMPORTANT -- provenance of the AM2 kinetic constants
%   ----------------------------------------------------
%   The AM2 values below are the published identification results of
%   Bernard et al. (2001) for an up-flow anaerobic fixed-bed reactor
%   treating raw industrial wine-distillery vinasse. They are NOT
%   substrate-specific to tomato residues co-digested with cattle manure.
%   They are used here as physically consistent starting values; the
%   methanogenic yield coefficient k6 (and optionally mu2max) should be
%   re-identified against measured data, or -- absent measurements --
%   rescaled so that the annual methane yield reproduces the site BMP used
%   by the static model. CALIBRATE_AM2_TO_BMP performs that rescaling and
%   must be reported as such in the manuscript.
%
%   References
%     Bernard, O., Hadj-Sadok, Z., Dochain, D., Genovesi, A., Steyer, J.-P.
%       (2001). Dynamical model development and parameter identification for
%       an anaerobic wastewater treatment process. Biotechnology and
%       Bioengineering, 75(4), 424-438.
%     Rosso, L., Lobry, J.R., Flandrois, J.P. (1993). An unexpected
%       correlation between cardinal temperatures of microbial growth
%       highlighted by a new model. J. Theor. Biol., 162(4), 447-463.
%     Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998). Crop
%       evapotranspiration. FAO Irrigation and Drainage Paper 56.
%     Monteith, J.L. (1977). Climate and the efficiency of crop production
%       in Britain. Phil. Trans. R. Soc. Lond. B, 281, 277-294.
%     Faiman, D. (2008). Assessing the outdoor operating temperature of
%       photovoltaic modules. Prog. Photovolt: Res. Appl., 16(4), 307-315.
%
%   See also APVAD_PARAMS, AM2_KINETICS, CALIBRATE_AM2_TO_BMP.

% =====================================================================
% 1. AM2 anaerobic digestion kinetics (Bernard et al. 2001)
% =====================================================================
D.am2.k1     = 42.14;   % g COD / g VSS      substrate degraded per acidogen
D.am2.k2     = 116.5;   % mmol VFA / g VSS   VFA produced by acidogenesis
D.am2.k3     = 268.0;   % mmol VFA / g VSS   VFA consumed by methanogenesis
D.am2.k4     = 50.6;    % mmol CO2 / g VSS   CO2 from acidogenesis
D.am2.k5     = 343.6;   % mmol CO2 / g VSS   CO2 from methanogenesis
D.am2.k6     = 453.0;   % mmol CH4 / g VSS   CH4 from methanogenesis

D.am2.mu1max = 1.20;    % 1/d   max acidogenic growth rate at T_ref
D.am2.KS1    = 7.10;    % g COD/L   Monod half-saturation (acidogenesis)
D.am2.mu2max = 0.74;    % 1/d   max methanogenic growth rate at T_ref
D.am2.KS2    = 9.28;    % mmol/L    Haldane half-saturation
D.am2.KI2    = 256.0;   % mmol/L    Haldane inhibition constant (= 16^2)

D.am2.alpha  = 0.50;    % -     fraction of biomass leaving with the effluent
                        %       (0 = perfect retention, 1 = ideal CSTR)
D.am2.kd1    = 0.02;    % 1/d   acidogen decay rate
D.am2.kd2    = 0.02;    % 1/d   methanogen decay rate

% Liquid-gas transfer and acid-base equilibrium
D.am2.kLa    = 19.8;    % 1/d   CO2 volumetric mass-transfer coefficient
D.am2.KH     = 16.0;    % mmol/(L atm)   Henry constant for CO2
D.am2.Pt     = 1.0;     % atm   total headspace pressure
D.am2.Kb     = 6.5e-7;  % mol/L bicarbonate equilibrium constant (for pH)

% Numerical floor: concentrations are clipped at this value inside the
% kinetics to keep the right-hand side smooth for the stiff solver.
D.am2.c_floor = 1e-9;

% =====================================================================
% 2. Temperature dependence of microbial growth
% =====================================================================
% 'arrhenius' reproduces the exponential form used by the static Python
% model (PARAMS.theta_arrhenius), capped at T_ref, so that the dynamic and
% static BMP responses are directly comparable. 'ctmi' is the cardinal
% temperature model with inflection (Rosso et al. 1993), which also
% captures inhibition above the optimum -- the physically correct choice
% when the digester can overheat (relevant at Ouagadougou).
D.temp.model    = 'arrhenius';   % 'arrhenius' | 'ctmi'
D.temp.T_min    = 5.0;           % degC   cardinal minimum (mesophilic)
D.temp.T_opt    = 37.0;          % degC   cardinal optimum
D.temp.T_max    = 45.0;          % degC   cardinal maximum
D.temp.T_ref    = 37.0;          % degC   reference for the Arrhenius form

% =====================================================================
% 3. Feed characterisation (tomato residue + cattle manure co-digestion)
% =====================================================================
D.feed.COD_per_VS   = 1.45;   % g COD / g VS   organic-matter conversion
D.feed.f_S1         = 0.70;   % -   biodegradable fraction entering S1
D.feed.S2in_mmolL   = 5.0;    % mmol/L   VFA already present in the feed
D.feed.Zin_mmolL    = 50.0;   % mmol/L   feed alkalinity (manure buffered)
D.feed.Cin_mmolL    = 55.0;   % mmol/L   feed total inorganic carbon
% Zin and Cin together fix the operating pH: at steady state the model
% gives pH = -log10(Kb*(C + S2 - Z)/(Z - S2)), so it is the DIFFERENCE
% Cin - Zin that matters, not either value alone. The pair above places the
% healthy operating point near pH 7.2, in the middle of the 7.0-7.8 band a
% stable mesophilic digester holds. Setting Cin = Zin, which looks harmless,
% drives the dissolved CO2 term to zero and the predicted pH above 8.
D.feed.T_in_degC    = 15.0;   % degC     feed temperature (= PARAMS.T_inlet)

% Residue handling: harvest delivers the season's residues over a short
% window; manure arrives continuously through the year.
D.feed.harvest_window_d = 10.0;  % d   spread of the residue delivery pulse

% Strength of the harvest removal, as the time-integral of the removal rate
% constant over the window. Removal is first-order in standing biomass, so
% a fraction exp(-strength) survives: a unit-area pulse would leave 37 % of
% the crop standing in the field, which is not a harvest. At 6 the residue
% left behind is 0.25 %.
D.feed.harvest_strength = 6.0;   % -   integral of the removal rate constant
D.feed.store_loss_1_d   = 0.001; % 1/d storage losses of stored VS
D.feed.min_store_kgVS   = 1.0;   % kg  store level below which feeding stops

% Target working buffer, expressed in days of feed. A farm digesting a
% once-a-year harvest necessarily holds several months of feedstock; four
% is a reasonable operating choice between harvests.
%
% This parameter does more work than it looks. The store obeys
%     dM/dt = supply - feed - k*M
% whose equilibrium is M* = (supply - feed)/k. With k = 0.001 1/d that is a
% time constant of a THOUSAND days -- by two orders of magnitude the
% slowest mode in the model. A store started away from M* is therefore
% still drifting after a decade of simulated years, and every "annual"
% methane figure along the way is partly fed by a draining inventory.
%
% Rather than integrate for decades or tune the feed rate iteratively, the
% equilibrium is imposed directly. Holding M* = buffer_days * feed/365
% makes annual storage losses exactly k*buffer_days*feed, so
%     feed = supply / (1 + k*buffer_days)
% and APVAD_INITIAL_STATE starts the store at the matching level. The store
% is then stationary from the first day and the slow mode never appears.
D.feed.buffer_days      = 120.0; % d   working stock held between harvests

% =====================================================================
% 4. Digester thermal dynamics
% =====================================================================
D.thermal.T_set_degC   = 37.0;   % degC   mesophilic setpoint
D.thermal.T_db_degC    = 1.0;    % degC   width of the controller dead band
D.thermal.use_biogas_boiler = true;  % allow biogas backup heating
D.thermal.boiler_eta   = 0.85;   % -      biogas boiler thermal efficiency
D.thermal.insulation_scale = 1.0;% -      multiplier on PARAMS.U_eff_W_m2K

% Heater capacity is SIZED, not fixed. The design heat load of a farm-scale
% digester is the wall loss at the coldest monthly mean plus the sensible
% heat to warm the influent -- a few hundred watts for a 4 m3 vessel, not
% the kilowatts a round number would suggest. BUILD_FORCING computes that
% load per design point and sets
%     Q_max = Q_max_factor * Q_design,   Kp = Q_max / prop_band_K.
% Getting this wrong is not a detail: an oversized heater holds the
% setpoint effortlessly in every scenario and destroys the model's ability
% to say anything about whether diverting f_PV_heat of the array is
% sufficient -- which is the question the coupling exists to answer.
D.thermal.Q_max_factor = 1.5;    % -      oversizing margin on design load
D.thermal.prop_band_K  = 3.0;    % K      error giving full heater output
D.thermal.Q_max_W      = [];     % W      set to a number to override sizing
D.thermal.Kp_W_per_K   = [];     % W/K    set to a number to override sizing

% =====================================================================
% 5. PV module thermal dynamics
% =====================================================================
% Lumped areal heat capacity of a glass/backsheet module. Typical values
% are 11-25 kJ/(m2 K), giving a thermal time constant of roughly 5-10 min.
D.pv.C_areal_kJ_m2K = 18.0;  % kJ/(m2 K)
D.pv.include_Pelec  = false; % subtract exported electrical power from the
                             % module energy balance (breaks the exact
                             % Faiman equivalence, so off by default)

% false -> algebraic Faiman, exactly the static model's expression
% true  -> first-order lag whose steady state IS that expression
%
% Default false, and the reason is a measured trade-off rather than a
% preference. The module time constant is about 12 minutes, two orders of
% magnitude faster than the digester's, so resolving it forces the stiff
% solver down to roughly one-minute steps for the whole year -- a full-year
% run goes from tens of seconds to several minutes. EXPERIMENT_PV_THERMAL_LAG
% quantifies what that buys: against hourly forcing the lag changes annual
% PV yield by a fraction of a percent, because the module spends the year
% alternating between lagging a rising sun and lagging a falling one and
% the two very nearly cancel.
%
% Turn it on for sub-hourly work, ramp-rate studies, or any question where
% the module's thermal memory is the subject rather than a detail.
D.pv.dynamic        = false;

% =====================================================================
% 6. Crop growth (Monteith RUE, thermal-time driven canopy)
% =====================================================================
D.crop.T_base_degC   = 10.0;   % degC   base temperature for tomato GDD
D.crop.T_cut_degC    = 30.0;   % degC   upper cut-off for GDD accumulation
D.crop.GDD_half      = 450.0;  % degC d thermal time at half canopy
D.crop.GDD_shape     = 3.0;    % -      logistic steepness of canopy build-up
D.crop.GDD_sen       = 1200.0; % degC d thermal time at onset of senescence
D.crop.LAI_max       = 3.5;    % m2/m2  maximum leaf area index
D.crop.k_ext         = 0.70;   % -      canopy light extinction coefficient
D.crop.RUE_g_MJ      = 2.50;   % g DM/MJ PAR  radiation use efficiency, an
                               % initial value only: CALIBRATE_CROP_RUE
                               % rescales it exactly against the site's
                               % reported open-field yield.
D.crop.DM_fruit      = 0.055;  % -      dry-matter content of fresh fruit

% Harvest index, derived rather than chosen, so that the crop model's
% residue output reconciles with the static model's residue accounting.
% The static model takes residue fresh mass as R_res = 0.80 times fruit
% fresh mass, at DM = 0.06 against the fruit's 0.055, so
%     residue_DM / fruit_DM = 0.80*0.06/0.055 = 0.873
%     HI = 1/(1 + 0.873)    = 0.534
% Using a round 0.6 instead would quietly starve the digester of about a
% fifth of its feedstock relative to the model being compared against.
D.crop.HI            = 0.534;  % -      fruit DM / total above-ground DM
D.crop.T_opt_growth  = 25.0;   % degC   optimum temperature for growth
D.crop.T_width_growth= 12.0;   % degC   half-width of the thermal response

% Between-season reset. Thermal time and standing biomass are states, and
% in a multi-year run they must return to zero before the next season or
% the second season starts on a canopy that is already senesced and
% produces almost nothing. Out of season, GDD relaxes to zero and any
% biomass still standing is cleared to the feedstock store, both with short
% time constants relative to the two or more months every site spends out
% of season.
D.crop.tau_gdd_reset_d = 2.0;  % d   thermal-time reset time constant
D.crop.tau_clear_d     = 5.0;  % d   clearance of residual standing biomass

% =====================================================================
% 7. Soil water balance (FAO-56 single coefficient)
% =====================================================================
D.water.Zr_m         = 0.80;   % m      effective rooting depth (tomato)
D.water.theta_fc     = 0.28;   % m3/m3  field capacity
D.water.theta_wp     = 0.13;   % m3/m3  permanent wilting point
D.water.p_depletion  = 0.40;   % -      readily available water fraction
D.water.Kc_ini       = 0.60;   % -      FAO-56 crop coefficient, initial
D.water.Kc_mid       = 1.15;   % -      mid-season
D.water.Kc_end       = 0.80;   % -      late season
D.water.irrigate     = true;   % apply automatic irrigation
D.water.irr_rate_mm_d= 12.0;   % mm/d   maximum irrigation application rate
D.water.rain_column  = '';     % name of a precipitation column if present
                               % in the climate CSV (empty -> assume 0)

% =====================================================================
% 8. Numerics
% =====================================================================
D.num.switch_width   = 0.25;   % smoothing width of the tanh switches used
                               % in place of hard if-statements, so that the
                               % right-hand side stays differentiable
D.num.rel_tol        = 1e-5;
D.num.abs_tol        = 1e-7;
D.num.max_step_d     = 1/24;   % d   never step over a forcing sample
D.num.nonnegative    = false;  % ODE15S NonNegative projection. The model
                               % already floors every concentration inside
                               % the kinetics, and the projection costs
                               % roughly a factor of two in runtime, so it
                               % is off by default; turn it on if a run ever
                               % produces a negative state.
end
