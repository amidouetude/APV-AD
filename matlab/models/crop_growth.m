function c = crop_growth(GDD, B_dm, PAR_AV, PAR_sat, T_day, inSeason, ...
                         Ks, harvestPulse, D)
%CROP_GROWTH Thermal-time canopy development and RUE biomass accumulation.
%
%   C = CROP_GROWTH(GDD, B_DM, PAR_AV, PAR_SAT, T_DAY, INSEASON, KS,
%   HARVESTPULSE, D) returns
%
%     C.dGDD_dt    thermal time derivative (degC d / d)
%     C.dB_dt      above-ground dry biomass derivative (kg DM/(ha d))
%     C.LAI        leaf area index (m2/m2)
%     C.f_i        fraction of PAR intercepted by the canopy (-)
%     C.APAR_MJ    absorbed PAR (MJ/(m2 d))
%     C.fT         thermal growth limitation (-)
%     C.harvest_kgDM_d  biomass leaving the field as residue (kg DM/(ha d))
%
%   Thermal time (growing degree days), accumulated in season only:
%       dGDD/dt = max(0, min(T_day, T_cut) - T_base)
%
%   Canopy development, logistic in thermal time with senescence:
%       LAI = LAI_max * GDD^n/(GDD_half^n + GDD^n) * 1/(1+exp((GDD-GDD_sen)/80))
%       f_i = 1 - exp(-k_ext * LAI)                        (Beer-Lambert)
%
%   Biomass accumulation (Monteith 1977 radiation use efficiency):
%       dB/dt = RUE * f_T * Ks * f_i * min(PAR_AV, PAR_sat) - harvest
%
%   Two things deserve comment.
%
%   First, min(PAR_AV, PAR_sat) is the same light-saturation cap the static
%   model applies, and it is what makes agrivoltaics work: above PAR_sat
%   the crop cannot use the extra light, so shading it costs nothing
%   agronomically while producing electricity. Applying the cap
%   instantaneously rather than to a season integral matters, because the
%   cap binds at midday and not at dawn.
%
%   Second, biomass is strictly linear in RUE: LAI depends on thermal time
%   alone, Ks on the water state alone, and the harvest term is linear in
%   B. So a single reference simulation is enough to calibrate RUE exactly
%   against the site's reported open-field yield -- no iteration. See
%   CALIBRATE_CROP_RUE.
%
%   References
%     Monteith, J.L. (1977). Phil. Trans. R. Soc. Lond. B, 281, 277-294.
%     McCree, K.J. (1971). Agricultural Meteorology, 9, 191-216.
%
%   See also SOIL_WATER, CALIBRATE_CROP_RUE.

cr = D.crop;

% ---- Thermal time -----------------------------------------------------
% Accumulates in season; relaxes back to zero out of season so that the
% next season starts on a bare field. Without the reset a multi-year run
% carries thermal time forward, the canopy begins the second season already
% past senescence, and the crop silently produces nothing -- while every
% RATIO in the results still looks plausible, because the reference run
% does the same thing.
dGDD_dt = max(min(T_day, cr.T_cut_degC) - cr.T_base_degC, 0) * inSeason ...
          - (1 - inSeason) * GDD / cr.tau_gdd_reset_d;

% ---- Canopy -----------------------------------------------------------
g   = max(GDD, 0);
build = g^cr.GDD_shape / (cr.GDD_half^cr.GDD_shape + g^cr.GDD_shape);
senesc = 1 / (1 + exp((g - cr.GDD_sen)/80));
LAI = cr.LAI_max * build * senesc * inSeason;
f_i = 1 - exp(-cr.k_ext * LAI);

% ---- Absorbed PAR -----------------------------------------------------
% The light-saturation cap is applied with a smoothed minimum. The crossing
% happens twice a day for the whole season, and a hard MIN puts a kink in
% the right-hand side at each one; the smoothing width, two percent of
% PAR_sat, is far below the precision the saturation point itself is known
% to (McCree's value for tomato is a fitted asymptote, not a sharp edge).
PAR_used = smooth_min(PAR_AV, PAR_sat, 0.02*PAR_sat);   % W/m2
PAR_used = max(PAR_used, 0);
APAR_W   = f_i * PAR_used;                    % W/m2
APAR_MJ  = APAR_W * 0.0864;                   % W/m2 -> MJ/(m2 d)

% ---- Thermal limitation on growth -------------------------------------
fT = exp(-((T_day - cr.T_opt_growth)/cr.T_width_growth)^2);

% ---- Biomass ----------------------------------------------------------
% RUE [g DM/MJ] * APAR [MJ/(m2 d)] = g/(m2 d); * 10 -> kg/(ha d)
growth = cr.RUE_g_MJ * APAR_MJ * fT * Ks * 10 * inSeason;

% Removal: the harvest itself, plus clearance of anything still standing
% once the season has ended. Both are first-order in standing biomass and
% both route to the feedstock store, so no dry matter is lost to
% bookkeeping.
removal_rate = harvestPulse + (1 - inSeason)/cr.tau_clear_d;   % 1/d
harvest = removal_rate * max(B_dm, 0);        % kg DM/(ha d)

c.dGDD_dt = dGDD_dt;
c.dB_dt   = growth - harvest;
c.LAI     = LAI;
c.f_i     = f_i;
c.APAR_MJ = APAR_MJ;
c.fT      = fT;
c.growth_kgDM_d  = growth;
c.harvest_kgDM_d = harvest;
end
