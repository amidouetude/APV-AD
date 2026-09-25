function C = compare_with_static(R_apv, R_ref, P, D, jsonPath)
%COMPARE_WITH_STATIC Benchmark the dynamic model against the Python results.
%
%   C = COMPARE_WITH_STATIC(R_APV, R_REF, P, D) compares the dynamic
%   simulation against the corresponding entry of
%   outputs/csv/optimal_params.json, the optimum found by the existing
%   Differential Evolution pipeline.
%
%   C = COMPARE_WITH_STATIC(..., JSONPATH) reads an explicit file.
%
%   Returned struct C holds a comparison table (C.table), the dynamic eLER
%   decomposition (C.dynamic) and the static reference (C.static).
%
%   How the eLER components are formed here
%   ---------------------------------------
%   Each component is the ratio of the APV run to the OPEN-FIELD REFERENCE
%   run of the same model, which is the same construction the static model
%   uses -- but every term is now an outcome rather than a ratio of
%   integrals:
%
%     LER_crop   = harvested dry matter (APV) / harvested dry matter (open)
%                  -- reflecting intercepted light, temperature and water
%                     stress through the season, not a PAR integral ratio
%     LER_PV     = PV exported (APV) / PV generated (open-field dense pack)
%     LER_biogas = methane energy (APV) / methane energy (open field)
%
%   Interpreting a disagreement
%   ---------------------------
%   Agreement on annual PV energy is expected and is a validation of the
%   port: the PV chain is the same physics, and with D.pv.dynamic = false
%   it is the same equations. Any material difference there means a bug,
%   not a discovery.
%
%   Disagreement on the biogas and crop terms is expected and is the
%   result. The static model converts an annual volatile-solids total
%   through a monthly-mean-temperature BMP and cannot represent a digester
%   that is cold in February, under-fed in June, or inhibited after the
%   harvest pulse. Report C.table alongside the reliability metrics in
%   R_APV.ANNUAL -- heat coverage, time below the mesophilic band, time
%   feed-limited -- since those are what explain the gap.
%
%   See also RUN_DYNAMIC_SIM, POSTPROCESS_SIM.

if nargin < 5 || isempty(jsonPath)
    here = fileparts(fileparts(mfilename('fullpath')));   % matlab/
    root = fileparts(here);
    jsonPath = fullfile(root, 'outputs', 'csv', 'optimal_params.json');
end

siteName = R_apv.site.name;

if isfile(jsonPath)
    raw = jsondecode(fileread(jsonPath));
    if isfield(raw, siteName)
        S = raw.(siteName);
    else
        warning('compare_with_static:siteMissing', ...
            'No entry for %s in %s.', siteName, jsonPath);
        S = struct();
    end
else
    warning('compare_with_static:jsonMissing', ...
        ['optimal_params.json not found at %s. Run 03_optimizer.ipynb, ' ...
         'or pass an explicit path. Only the dynamic results are ' ...
         'reported.'], jsonPath);
    S = struct();
end

A = R_apv.annual;
B = R_ref.annual;

% ---------------------------------------------------------------------
% Dynamic eLER decomposition
% ---------------------------------------------------------------------
dyn.LER_crop   = A.harvest_DM_kg_ha / max(B.harvest_DM_kg_ha, eps);
dyn.LER_PV     = A.PV_sold_kWh_ha   / max(B.PV_total_kWh_ha,  eps);
dyn.LER_biogas = A.biogas_kWh_ha    / max(B.biogas_kWh_ha,    eps);
dyn.eLER       = dyn.LER_crop + dyn.LER_PV + dyn.LER_biogas;
dyn.LER_2C     = dyn.LER_crop + dyn.LER_PV;

dyn.PV_total_MWh_ha    = A.PV_total_MWh_ha;
dyn.PV_sold_MWh_ha     = A.PV_sold_MWh_ha;
dyn.biogas_total_MWh_ha= A.biogas_MWh_ha;
dyn.heat_demand_MWh_ha = A.heat_demand_kWh_ha/1000;
dyn.BMP_eff_NmL_gVS    = A.BMP_eff_NmL_gVS;
dyn.OLR_kgVS_m3d       = A.OLR_mean_kgVS_m3d;
dyn.T_dig_avg_degC     = A.T_dig_mean_degC;
dyn.GCR_pct            = 100*A.GCR;
dyn.ESR                = A.ESR;
dyn.LER_water_season_pct = A.LER_water_season_pct;
dyn.ET0_season_mm      = A.ET0_season_mm;
dyn.W_saved_mm_season  = A.W_saved_static_mm;

% ---------------------------------------------------------------------
% Side-by-side table
% ---------------------------------------------------------------------
dyn.LER_crop_PARstyle = A.LER_crop_PARstyle;

pairs = {
  'PV_total_MWh_ha',      'PV_total_MWh_ha',      dyn.PV_total_MWh_ha
  'PV_sold_MWh_ha',       'PV_sold_MWh_ha',       dyn.PV_sold_MWh_ha
  'biogas_total_MWh_ha',  'biogas_total_MWh_ha',  dyn.biogas_total_MWh_ha
  'heat_demand_MWh_ha',   'heat_demand_MWh_ha',   dyn.heat_demand_MWh_ha
  'BMP_avg_NmL_gVS',      'BMP_eff_NmL_gVS',      dyn.BMP_eff_NmL_gVS
  'OLR_kgVS_m3d',         'OLR_kgVS_m3d',         dyn.OLR_kgVS_m3d
  'T_dig_avg_degC',       'T_dig_avg_degC',       dyn.T_dig_avg_degC
  'GCR_pct',              'GCR_pct',              dyn.GCR_pct
  'LER_crop',             'LER_crop (PAR-ratio defn)', dyn.LER_crop_PARstyle
  'LER_crop',             'LER_crop (biomass defn)',   dyn.LER_crop
  'LER_PV_defA',          'LER_PV',               dyn.LER_PV
  'LER_biogas',           'LER_biogas',           dyn.LER_biogas
  'eLER_defA',            'eLER',                 dyn.eLER
  'ET0_season_mm',        'ET0_season_mm',        dyn.ET0_season_mm
  'W_saved_mm_season',    'W_saved_mm_season',    dyn.W_saved_mm_season
  'LER_water_season_pct', 'LER_water_season_pct', dyn.LER_water_season_pct
};

n = size(pairs,1);
Quantity = strings(n,1);
Static   = nan(n,1);
Dynamic  = nan(n,1);

for i = 1:n
    Quantity(i) = string(pairs{i,2});
    Dynamic(i)  = pairs{i,3};
    if isfield(S, pairs{i,1})
        Static(i) = S.(pairs{i,1});
    end
end

Diff_pct = 100*(Dynamic - Static)./abs(Static);
C.table = table(Quantity, Static, Dynamic, Diff_pct);

C.dynamic = dyn;
C.static  = S;
C.site    = siteName;
C.reliability = struct( ...
    'heat_coverage',            A.heat_coverage, ...
    'pv_share_of_heat',         A.pv_share_of_heat, ...
    'frac_time_below_T_min',    A.frac_time_below_min, ...
    'frac_time_mesophilic',     A.frac_time_mesophilic, ...
    'frac_time_feed_limited',   A.frac_time_feed_limited, ...
    'frac_time_VFA_inhibited',  A.frac_time_inhibited, ...
    'frac_season_water_stressed', A.frac_season_stressed, ...
    'pH_min',                   A.pH_min, ...
    'x_CH4_mean',               A.x_CH4_mean);

% ---------------------------------------------------------------------
% Report
% ---------------------------------------------------------------------
fprintf('\n=== %s: dynamic (MATLAB) vs static (Python DE optimum) ===\n', ...
        siteName);
disp(C.table);

fprintf('--- Reliability metrics available only from the dynamic model ---\n');
fprintf('  heat demand covered            : %6.1f %%\n', 100*A.heat_coverage);
fprintf('  of delivered heat, from PV     : %6.1f %%\n', 100*A.pv_share_of_heat);
fprintf('  time below %4.1f degC           : %6.1f %%\n', ...
        P.T_dig_min_degC, 100*A.frac_time_below_min);
fprintf('  time in mesophilic band (>=35) : %6.1f %%\n', ...
        100*A.frac_time_mesophilic);
fprintf('  time feed-limited              : %6.1f %%\n', ...
        100*A.frac_time_feed_limited);
fprintf('  time VFA above inhibition      : %6.1f %%  (S2 > %.0f mmol/L)\n', ...
        100*A.frac_time_inhibited, A.S2_crit_mmolL);
fprintf('  season under water stress      : %6.1f %%\n', ...
        100*A.frac_season_stressed);
fprintf('  minimum pH                     : %6.2f\n', A.pH_min);
fprintf('  mean methane fraction          : %6.1f %% (static assumes %.0f %%)\n', ...
        100*A.x_CH4_mean, 100*P.CH4_fraction);
fprintf('  digester thermal time constant : %6.2f d\n', A.tau_dig_d);

fprintf('\n--- Periodicity of the reported year ---\n');
pc = A.periodicity;
fprintf('  worst state drift              : %6.2f %% (%s)\n', ...
        100*pc.worst_rel, pc.worst_state);
fprintf('  feedstock store drift          : %+6.0f kg VS over the year\n', ...
        pc.store_drift_kgVS);
if pc.converged
    fprintf('  -> periodic steady state reached\n');
else
    fprintf(['  -> NOT at periodic steady state. Part of the reported ' ...
             'annual output came from\n     draining an inventory rather ' ...
             'than from that year''s feedstock. Increase\n     ' ...
             'opts.n_years before quoting these figures.\n']);
end
fprintf('\n');
end
