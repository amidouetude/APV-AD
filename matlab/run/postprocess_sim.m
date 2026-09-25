function R = postprocess_sim(R, D, P)
%POSTPROCESS_SIM Recover algebraic outputs and annual aggregates.
%
%   R = POSTPROCESS_SIM(R, D, P) evaluates APVAD_RHS a second time on the
%   output grid, this time requesting its algebraic outputs, and fills
%
%     R.out      struct of n x 1 time series (powers, flows, pH, stresses)
%     R.annual   struct of annual totals, averages and reliability metrics
%
%   Annual energy, water and biomass totals are taken as DIFFERENCES OF THE
%   ACCUMULATOR STATES across the reported year, not by re-integrating the
%   sampled output series. The accumulators were integrated by the solver
%   under its own error control at its own step size; re-integrating an
%   hourly resample would silently add a quadrature error of its own, and
%   it is precisely the sharp signals -- a heater switching, a harvest
%   pulse -- that an hourly trapezoid handles worst.
%
%   See also RUN_DYNAMIC_SIM, APVAD_RHS, COMPARE_WITH_STATIC.

ix = apvad_states();
F  = R.F;
t  = R.t_abs;
n  = numel(t);

% ---------------------------------------------------------------------
% 1. Algebraic outputs on the output grid
% ---------------------------------------------------------------------
[~, o1] = apvad_rhs(t(1), R.X(1,:)', F, D, P);
flds = fieldnames(o1);
for k = 1:numel(flds)
    out.(flds{k}) = zeros(n, 1);
    out.(flds{k})(1) = o1.(flds{k});
end
for i = 2:n
    [~, oi] = apvad_rhs(t(i), R.X(i,:)', F, D, P);
    for k = 1:numel(flds)
        out.(flds{k})(i) = oi.(flds{k});
    end
end

% States exposed alongside the algebraic outputs, for convenience.
for k = 1:ix.n
    out.(ix.names{k}) = R.X(:, k);
end
R.out = out;

% ---------------------------------------------------------------------
% 2. Annual totals from the accumulators
% ---------------------------------------------------------------------
dAcc = @(j) R.X(end, j) - R.X(1, j);

A.PV_total_kWh_ha   = dAcc(ix.E_pv);
A.PV_sold_kWh_ha    = dAcc(ix.E_pv_sold);
A.PV_heat_kWh_ha    = A.PV_total_kWh_ha - A.PV_sold_kWh_ha;
A.biogas_kWh_ha     = dAcc(ix.E_ch4);
A.heat_delivered_kWh_ha = dAcc(ix.E_heat);
A.biogas_to_boiler_kWh_ha = dAcc(ix.E_ch4_heat);
A.irrigation_mm     = dAcc(ix.W_irr);
A.harvest_DM_kg_ha  = dAcc(ix.M_harv);

A.PV_total_MWh_ha   = A.PV_total_kWh_ha/1000;
A.PV_sold_MWh_ha    = A.PV_sold_kWh_ha/1000;
A.biogas_MWh_ha     = A.biogas_kWh_ha/1000;

% Net biogas available after the digester has heated itself.
A.biogas_net_kWh_ha = A.biogas_kWh_ha - A.biogas_to_boiler_kWh_ha;

% ---------------------------------------------------------------------
% 3. Crop and feedstock
% ---------------------------------------------------------------------
A.fruit_DM_kg_ha    = D.crop.HI * A.harvest_DM_kg_ha;
A.fruit_fresh_t_ha  = A.fruit_DM_kg_ha / D.crop.DM_fruit / 1000;
A.residue_DM_kg_ha  = (1 - D.crop.HI) * A.harvest_DM_kg_ha;
A.residue_VS_kg_ha  = A.residue_DM_kg_ha * R.site.VS_fraction;
A.LAI_max           = max(out.LAI);
A.APAR_season_MJ    = trapz(t, out.APAR_MJ);

% ---------------------------------------------------------------------
% 4. Digester operation and health
% ---------------------------------------------------------------------
season = out.in_season > 0.5;

% The static model's own LER_crop numerator, recomputed here on the dynamic
% run so that the two definitions can be compared directly rather than
% confounded with a modelling difference. The static formula is
%     LER_crop = sum(min(PAR_AV, PAR_sat)) / sum(PAR_open)
% over the growing season. Note that the numerator is capped at light
% saturation while the denominator is not. That asymmetry means the ratio
% falls below one even with NO shading at all, so it charges the
% agrivoltaic system for the crop's own light saturation -- an effect that
% has nothing to do with the panels. COMPARE_WITH_STATIC reports this
% quantity beside the biomass-based ratio precisely so the definitional
% difference is not mistaken for a physical one.
PAR_open_dyn = out.PAR_AV ./ max(1 - out.F_shad, 1e-9);
A.LER_crop_PARstyle = trapz(t(season), ...
                        min(out.PAR_AV(season), R.site.PAR_sat)) ...
                      / max(trapz(t(season), PAR_open_dyn(season)), eps);

A.T_dig_mean_degC   = mean(out.T_dig);
A.T_dig_min_degC    = min(out.T_dig);
A.T_dig_max_degC    = max(out.T_dig);
A.frac_time_mesophilic = mean(out.T_dig >= 35);
A.frac_time_below_min  = mean(out.T_dig <  P.T_dig_min_degC);

A.pH_mean           = mean(out.pH);
A.pH_min            = min(out.pH);
A.frac_time_pH_low  = mean(out.pH < 6.8);

A.VFA_mean_mmolL    = mean(out.S2);
A.VFA_max_mmolL     = max(out.S2);
% Haldane inhibition sets in above sqrt(KS2*KI2); time spent past it is a
% direct measure of how close the reactor runs to souring.
S2_crit = sqrt(D.am2.KS2 * D.am2.KI2);
A.S2_crit_mmolL     = S2_crit;
A.frac_time_inhibited = mean(out.S2 > S2_crit);

A.OLR_mean_kgVS_m3d = mean(out.OLR_kgVS_m3d);
A.OLR_max_kgVS_m3d  = max(out.OLR_kgVS_m3d);
A.x_CH4_mean        = mean(out.x_CH4(out.Q_biogas_Nm3_d > 1e-6));
A.CH4_Nm3_ha        = trapz(t, out.Q_CH4_Nm3_d);
A.BMP_eff_NmL_gVS   = 1e6 * A.CH4_Nm3_ha ...
                      / max(trapz(t, out.VS_feed_kg_d)*1000, eps);
A.frac_time_feed_limited = mean(out.store_avail < 0.99);
A.store_min_kgVS    = min(out.M_res);
A.store_max_kgVS    = max(out.M_res);

% ---------------------------------------------------------------------
% 5. Heating reliability -- the metric the static model cannot produce
% ---------------------------------------------------------------------
A.heat_demand_kWh_ha  = trapz(t, out.Q_demand_W) * 24/1000;
A.heat_deficit_kWh_ha = trapz(t, out.deficit_W)  * 24/1000;
A.heat_from_pv_kWh_ha = trapz(t, out.Q_pv_W)     * 24/1000;
A.heat_from_boiler_kWh_ha = trapz(t, out.Q_boiler_W) * 24/1000;
A.heat_coverage       = 1 - A.heat_deficit_kWh_ha ...
                        / max(A.heat_demand_kWh_ha, eps);
A.pv_share_of_heat    = A.heat_from_pv_kWh_ha ...
                        / max(A.heat_delivered_kWh_ha, eps);
A.frac_time_heat_short = mean(out.deficit_W > 1);

% Energy self-sufficiency, defined as in the static model.
A.ESR = (A.PV_total_kWh_ha + A.biogas_kWh_ha) ...
        / max(A.heat_demand_kWh_ha, 1);

% ---------------------------------------------------------------------
% 6. Water
% ---------------------------------------------------------------------
A.ET0_annual_mm     = trapz(t, out.ET0);
A.ET0_season_mm     = trapz(t(season), out.ET0(season));
A.ETc_season_mm     = trapz(t(season), out.ETc_mm_d(season));
A.ETc_open_season_mm= trapz(t(season), out.ETc_open(season));
A.W_saved_dyn_mm    = A.ETc_open_season_mm - A.ETc_season_mm;
% The static model's definition, retained for a like-for-like comparison:
% W_saved = sum(alpha_shade * F_shad * ET0) over the season, with no crop
% coefficient and no water-stress feedback.
A.W_saved_static_mm = trapz(t(season), ...
                        P.alpha_shade * out.F_shad(season) .* out.ET0(season));
A.LER_water_season_pct = 100 * A.W_saved_static_mm / max(A.ET0_season_mm, eps);
A.Ks_mean_season    = mean(out.Ks(season));
A.frac_season_stressed = mean(out.Ks(season) < 0.99);
A.TAW_mm            = out.TAW(1);
A.RAW_mm            = out.RAW(1);
A.Dr_max_mm         = max(out.Dr);

% ---------------------------------------------------------------------
% 7. PV thermal
% ---------------------------------------------------------------------
lit = out.G_eff > 20;
A.T_cell_mean_lit_degC = mean(out.T_cell(lit));
A.T_cell_max_degC      = max(out.T_cell);
A.T_cell_lag_rms_K     = sqrt(mean((out.T_cell(lit) - out.T_cell_ss(lit)).^2));
A.GCR                  = F.GCR;
A.n_modules            = F.n_modules;

% ---------------------------------------------------------------------
% 8. Periodicity check
% ---------------------------------------------------------------------
% A reported year is only meaningful if the system returns to where it
% started -- otherwise part of the "annual" methane came out of a draining
% feedstock store or a shrinking biomass inventory rather than out of that
% year's feedstock, and the reported specific yield is inflated by exactly
% that drift. The check is on the non-accumulator states only, since the
% accumulators are meant to grow.
% Drift is judged against how far the state MOVES within the year, not
% against its starting value. A state that swings seasonally between zero
% and several tonnes -- standing biomass, the feedstock store -- would
% otherwise report an enormous relative drift for a physically trivial
% offset, or none at all if it happens to start near its own maximum.
dyn_states = setdiff(1:ix.n, ix.accumulators);
x_start = R.X(1,   dyn_states);
x_end   = R.X(end, dyn_states);
swing   = max(R.X(:, dyn_states), [], 1) - min(R.X(:, dyn_states), [], 1);
scale   = max([swing; abs(x_start); 1e-6*ones(1,numel(dyn_states))], [], 1);

A.periodicity.names    = ix.names(dyn_states);
A.periodicity.rel_drift = (x_end - x_start) ./ scale;
[A.periodicity.worst_rel, iw] = max(abs(A.periodicity.rel_drift));
A.periodicity.worst_state = ix.names{dyn_states(iw)};
A.periodicity.store_drift_kgVS = R.X(end, ix.M_res) - R.X(1, ix.M_res);
A.periodicity.converged = A.periodicity.worst_rel < 0.05;

A.tau_dig_d = out.tau_dig_d(1);
A.mode      = R.mode;
A.site      = R.site.name;

R.annual = A;
end
