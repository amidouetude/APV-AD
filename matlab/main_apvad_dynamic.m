function OUT = main_apvad_dynamic(siteName, opts)
%MAIN_APVAD_DYNAMIC Top-level driver for the dynamic APV-AD model.
%
%   OUT = MAIN_APVAD_DYNAMIC()            runs Konya with the DE optimum
%   OUT = MAIN_APVAD_DYNAMIC(SITE)        runs one of the four study sites
%   OUT = MAIN_APVAD_DYNAMIC(SITE, OPTS)  overrides defaults
%
%   OPTS fields
%     design       struct of the six design variables; default is the DE
%                  optimum for the site from outputs/csv/optimal_params.json
%     n_years      years integrated including spin-up (default 2)
%     calibrate    true (default) -- anchor RUE to the reported open-field
%                  yield and k6 to the reference BMP assay
%     pv_dynamic   override D.pv.dynamic
%     temp_model   'arrhenius' (default) or 'ctmi'
%     plots        true (default) -- draw the dynamics figure
%     save_results true (default) -- write outputs/csv/dynamic_<site>.mat
%                  and the summary CSV
%
%   OUT contains
%     OUT.apv       APV simulation (states, outputs, annual aggregates)
%     OUT.ref       open-field reference simulation
%     OUT.compare   comparison against the static Python optimum
%     OUT.D, OUT.P  the parameter sets actually used, after calibration
%     OUT.calib     calibration report
%
%   Pipeline
%   --------
%     1  load the site's hourly TMY climate
%     2  reference (open-field) run, one year, to calibrate RUE exactly
%     3  calibrate k6 against the BMP assay condition (closed form)
%     4  reference run with calibrated parameters -> LER denominators
%     5  APV run
%     6  compare, report, plot
%
%   Steps 2 and 4 are separate runs because rescaling RUE changes residue
%   production, which changes the digester feed, which changes everything
%   downstream. RUE itself is recovered exactly from a single run because
%   biomass is linear in it; the second run propagates the consequence.
%
%   Example
%     out = main_apvad_dynamic('Konya');
%     out = main_apvad_dynamic('Ouagadougou', struct('temp_model','ctmi'));
%
%   See also RUN_DYNAMIC_SIM, COMPARE_WITH_STATIC, PLOT_APVAD_DYNAMICS.

% ---------------------------------------------------------------------
% Path setup: make every subfolder of matlab/ visible.
% ---------------------------------------------------------------------
thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(thisDir));

if nargin < 1 || isempty(siteName), siteName = 'Konya'; end
if nargin < 2, opts = struct(); end
opts = setdef(opts, 'n_years',      2);
opts = setdef(opts, 'calibrate',    true);
opts = setdef(opts, 'temp_model',   'arrhenius');
opts = setdef(opts, 'plots',        true);
opts = setdef(opts, 'save_results', true);
opts = setdef(opts, 'verbose',      true);

projRoot = fileparts(thisDir);

% ---------------------------------------------------------------------
% Parameters and design point
% ---------------------------------------------------------------------
P    = apvad_params();
D    = apvad_dyn_params();
site = apvad_sites(siteName);

D.temp.model = opts.temp_model;
if isfield(opts, 'pv_dynamic') && ~isempty(opts.pv_dynamic)
    D.pv.dynamic = opts.pv_dynamic;
end

if isfield(opts, 'design') && ~isempty(opts.design)
    design = opts.design;
else
    design = load_optimal_design(siteName, projRoot);
end

fprintf('\n==================================================\n');
fprintf(' APV-AD dynamic model -- %s\n', siteName);
fprintf('==================================================\n');
fprintf(' design: beta %.2f deg | d_row %.2f m | H %.2f m\n', ...
        design.beta_deg, design.d_row_m, design.H_m_m);
fprintf('         V_dig %.2f m3 | HRT %.1f d | f_PV_heat %.3f\n', ...
        design.V_dig_m3, design.HRT_days, design.f_PV_heat);
fprintf(' growth temperature model: %s | PV thermal dynamic: %d\n\n', ...
        D.temp.model, D.pv.dynamic);

% ---------------------------------------------------------------------
% 1. Climate
% ---------------------------------------------------------------------
clim = load_climate(siteName, fullfile(projRoot, 'outputs', 'csv'));

% ---------------------------------------------------------------------
% 2-3. Calibration
% ---------------------------------------------------------------------
calib = struct();
if opts.calibrate
    fprintf('--- calibration ---------------------------------\n');

    % RUE against the reported open-field yield (exact, one run).
    R_cal = run_dynamic_sim(clim, site, design, P, D, ...
                struct('mode','reference','n_years',1,'verbose',opts.verbose));
    [D, calib.rue] = calibrate_crop_rue(R_cal, D);
    fprintf(['  RUE  %.3f -> %.3f g DM/MJ  (open-field yield %.1f -> ' ...
             '%.1f t/ha)\n'], calib.rue.RUE_before, calib.rue.RUE_after, ...
             calib.rue.yield_sim_t_ha, calib.rue.yield_target_t_ha);

    % k6 against the BMP assay condition (closed form, no run).
    Ftmp = build_forcing(clim, site, design, P, D, 'apv');
    [D, calib.bmp] = calibrate_am2_to_bmp(D, P, design.HRT_days, ...
                                          Ftmp.c_VS_in_gL);
    fprintf(['  k6   %.1f -> %.1f mmol CH4/g VSS  (BMP at 37 degC ' ...
             '%.1f -> %.1f NmL/g VS)\n\n'], calib.bmp.k6_before, ...
             calib.bmp.k6_after, calib.bmp.BMP_model, calib.bmp.BMP_after);
end

% ---------------------------------------------------------------------
% 4. Open-field reference run (LER denominators)
% ---------------------------------------------------------------------
% The reference case gets the same stationary-store treatment. Its residue
% supply is the site's reported open-field figure by construction, since
% that is exactly what RUE was just calibrated to reproduce.
VS_ref_supply = site.fruit_yield_t_ha*1000*site.R_res*site.DM*site.VS_fraction ...
              + site.manure_input_t_ha*1000*site.manure_DM*site.manure_VS_fraction;
design_ref = design;
design_ref.VS_ann_kg = VS_ref_supply ...
                       / (1 + D.feed.store_loss_1_d*D.feed.buffer_days);

R_ref = run_dynamic_sim(clim, site, design_ref, P, D, ...
            struct('mode','reference','n_years',opts.n_years, ...
                   'verbose',opts.verbose));

% ---------------------------------------------------------------------
% 5. Agrivoltaic run: converge the feedstock basis, then simulate
% ---------------------------------------------------------------------
% The digester is fed at a constant volumetric rate; what has to be chosen
% is the influent concentration, i.e. how much volatile solids the design
% assumes are available each year. Two things make the obvious choice
% wrong. Under an array the crop yields less residue than the open-field
% figure the site parameters report, and the feedstock store loses about
% 0.1 % of its contents per day to degradation.
%
% That loss coefficient gives the store a time constant of roughly a
% THOUSAND days. It is the slowest mode in the whole model by two orders of
% magnitude, and it means a store that starts out of balance is still out
% of balance after three simulated years -- the annual methane figure would
% then include feedstock drawn from an inventory rather than produced that
% year, and would drift downwards for a decade.
%
% So the basis is converged rather than guessed: run a year, measure how
% far the store drifted, and correct the basis by exactly that amount. The
% store balance dM/dt = supply - feed - k*M is linear in the feed rate, so
% this is a Newton step on a linear function and lands in two or three
% passes. Convergence means the reported year is genuinely fed by that
% year's harvest.
% Pass A: one year on the nominal open-field basis, purely to measure how
% much residue the crop actually delivers under the array. The crop
% sub-model does not depend on the digester at all, so a single pass
% determines the supply exactly -- no iteration is required.
R_apvA = run_dynamic_sim(clim, site, design, P, D, ...
            struct('mode','apv','n_years',1,'verbose',false));

supply = R_apvA.annual.residue_VS_kg_ha + R_apvA.F.manure_VS_kg_d*365;

% Feed rate consistent with a stationary store (see D.feed.buffer_days):
%     feed = supply / (1 + k*buffer_days)
% the remainder being lost to degradation in storage.
k_loss = D.feed.store_loss_1_d;
basis  = supply / (1 + k_loss*D.feed.buffer_days);

calib.VS_ann_nominal  = R_apvA.F.VS_ann_nom_kg;
calib.VS_ann_supply   = supply;
calib.VS_ann_basis    = basis;
calib.VS_storage_loss = supply - basis;

fprintf(['--- feedstock: nominal %.0f -> realised supply %.0f kg VS/ha ' ...
         '(%.0f %% of open field)\n'], calib.VS_ann_nominal, supply, ...
         100*supply/max(calib.VS_ann_nominal, eps));
fprintf(['--- feed basis %.0f kg VS/ha/yr, storage losses %.0f kg ' ...
         '(%.1f %%), buffer %.0f d\n\n'], basis, calib.VS_storage_loss, ...
         100*calib.VS_storage_loss/max(supply,eps), D.feed.buffer_days);

design_apv = design;
design_apv.VS_ann_kg = basis;

R_apv = run_dynamic_sim(clim, site, design_apv, P, D, ...
            struct('mode','apv','n_years',opts.n_years, ...
                   'verbose',opts.verbose));

% ---------------------------------------------------------------------
% 6. Compare and report
% ---------------------------------------------------------------------
C = compare_with_static(R_apv, R_ref, P, D, ...
        fullfile(projRoot,'outputs','csv','optimal_params.json'));

OUT.apv     = R_apv;
OUT.ref     = R_ref;
OUT.compare = C;
OUT.D       = D;
OUT.P       = P;
OUT.site    = site;
OUT.design  = design;
OUT.calib   = calib;

% ---------------------------------------------------------------------
% 7. Figures and files
% ---------------------------------------------------------------------
if opts.plots
    figDir = fullfile(projRoot, 'outputs', 'figures');
    if ~isfolder(figDir), mkdir(figDir); end
    OUT.fig = plot_apvad_dynamics(OUT, figDir);
end

if opts.save_results
    csvDir = fullfile(projRoot, 'outputs', 'csv');
    if ~isfolder(csvDir), mkdir(csvDir); end

    % Save the whole result struct, minus the raw solver output and the
    % forcing interpolants, so PLOT_FROM_SAVED can redraw every figure
    % later without repeating the integration.
    OUTs = OUT;
    OUTs.apv = strip_solution(OUT.apv);
    OUTs.ref = strip_solution(OUT.ref);

    matFile = fullfile(csvDir, sprintf('dynamic_%s.mat', siteName));
    save(matFile, 'OUTs', '-v7.3');

    csvFile = fullfile(csvDir, sprintf('dynamic_%s_summary.csv', siteName));
    writetable(C.table, csvFile);

    fprintf('saved: %s\n', matFile);
    fprintf('saved: %s\n', csvFile);
end
end

% =====================================================================
function s = setdef(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end

% =====================================================================
function design = load_optimal_design(siteName, projRoot)
%LOAD_OPTIMAL_DESIGN Read the DE optimum, or fall back to a sane default.
jsonPath = fullfile(projRoot, 'outputs', 'csv', 'optimal_params.json');
fields = {'beta_deg','d_row_m','H_m_m','V_dig_m3','HRT_days','f_PV_heat'};

if isfile(jsonPath)
    raw = jsondecode(fileread(jsonPath));
    if isfield(raw, siteName)
        S = raw.(siteName);
        if all(cellfun(@(f) isfield(S,f), fields))
            for k = 1:numel(fields)
                design.(fields{k}) = S.(fields{k});
            end
            return
        end
    end
end

warning('main_apvad_dynamic:noOptimum', ...
    ['No optimum for %s in optimal_params.json; using the midpoint of ' ...
     'the DE bounds. Run 03_optimizer.ipynb for the design point the ' ...
     'manuscript reports.'], siteName);
P = apvad_params();
mid = 0.5*(P.design_lb + P.design_ub);
for k = 1:numel(fields)
    design.(fields{k}) = mid(k);
end
end

% =====================================================================
function R = strip_solution(R)
%STRIP_SOLUTION Drop the raw solver structure and forcing interpolants
%   before saving. Both are large, neither is needed to read the results,
%   and function handles inside griddedInterpolant objects do not survive a
%   round trip cleanly.
if isfield(R, 'sol'), R = rmfield(R, 'sol'); end
if isfield(R, 'F')
    keep = {'mode','design','n_modules','GCR','d_row_m','A_pv_m2', ...
            'doy_plant','doy_harvest','season_len_d','V_dig_m3', ...
            'Q_in_m3_d','A_wall_m2','c_VS_in_gL','manure_VS_kg_d', ...
            'VS_ann_nom_kg','t_period'};
    Fs = struct();
    for k = 1:numel(keep)
        if isfield(R.F, keep{k}), Fs.(keep{k}) = R.F.(keep{k}); end
    end
    R.F = Fs;
end
end
