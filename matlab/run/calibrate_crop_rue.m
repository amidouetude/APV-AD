function [D, info] = calibrate_crop_rue(R_ref, D)
%CALIBRATE_CROP_RUE Scale radiation use efficiency to the reported yield.
%
%   [D, INFO] = CALIBRATE_CROP_RUE(R_REF, D) rescales D.crop.RUE_g_MJ so
%   that the OPEN-FIELD reference simulation R_REF reproduces the fresh
%   fruit yield recorded for the site in APVAD_SITES (60 t/ha at the three
%   warm sites, 30 t/ha at Freiburg).
%
%   The rescaling is exact and needs no iteration. Canopy development
%   depends on thermal time alone, the water-stress coefficient on the soil
%   water state alone, and the harvest term is linear in standing biomass,
%   so biomass is strictly proportional to RUE. One reference run therefore
%   determines the constant exactly:
%
%       RUE_new = RUE_old * yield_reported / yield_simulated
%
%   Anchoring the crop model to the same yield the static model assumes is
%   what makes LER_crop comparable between the two. The static model takes
%   the site yield as given and reduces it by a PAR integral ratio; here
%   the reduction emerges from intercepted light, temperature and water
%   stress acting through the season. If the two disagree, the disagreement
%   is about the mechanism, not about the baseline.
%
%   The reference run must be produced with mode 'reference' -- an
%   open-field, dense-packed array with no inter-row ground shading.
%
%   INFO reports the simulated and target yields and the scale factor.
%
%   See also RUN_DYNAMIC_SIM, CROP_GROWTH, BUILD_FORCING.

if ~strcmpi(R_ref.mode, 'reference')
    warning('calibrate_crop_rue:notReference', ...
        ['Calibrating against a run in mode "%s". RUE should be anchored ' ...
         'to the OPEN-FIELD reference yield; using a shaded run will ' ...
         'inflate RUE to compensate for the shading.'], R_ref.mode);
end

y_sim    = R_ref.annual.fruit_fresh_t_ha;
y_target = R_ref.site.fruit_yield_t_ha;

if ~(y_sim > 0)
    error('calibrate_crop_rue:zeroYield', ...
        ['The reference simulation produced no fruit yield, so RUE cannot ' ...
         'be scaled. Check the growing-season mask (in_season) and the ' ...
         'harvest date for %s.'], R_ref.site.name);
end

scale = y_target / y_sim;

info.RUE_before      = D.crop.RUE_g_MJ;
info.yield_sim_t_ha  = y_sim;
info.yield_target_t_ha = y_target;
info.scale           = scale;

D.crop.RUE_g_MJ = D.crop.RUE_g_MJ * scale;
info.RUE_after  = D.crop.RUE_g_MJ;
end
