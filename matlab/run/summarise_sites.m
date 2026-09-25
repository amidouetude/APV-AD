function T = summarise_sites(siteNames, csvDir)
%SUMMARISE_SITES Cross-site table built from saved dynamic-model results.
%
%   T = SUMMARISE_SITES() reads outputs/csv/dynamic_<site>.mat for the four
%   study sites, assembles the comparison table and writes
%   outputs/csv/dynamic_all_sites.csv.
%
%   T = SUMMARISE_SITES(SITENAMES, CSVDIR) for a subset or another folder.
%
%   Kept separate from RUN_ALL_SITES so the summary can be rebuilt, or its
%   columns changed, without repeating four multi-minute integrations. Sites
%   with no saved result are skipped with a warning rather than failing the
%   whole table.
%
%   The PERIODIC column is the one to read first: any site marked false did
%   not reach a periodic steady state, so part of its reported annual output
%   came from draining an inventory and its figures are not quotable.
%
%   See also RUN_ALL_SITES, MAIN_APVAD_DYNAMIC, POSTPROCESS_SIM.

thisDir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(thisDir));
projRoot = fileparts(thisDir);

if nargin < 1 || isempty(siteNames)
    siteNames = {'Konya','Almeria','Ouagadougou','Freiburg'};
end
if nargin < 2 || isempty(csvDir)
    csvDir = fullfile(projRoot,'outputs','csv');
end

rows = {};
for k = 1:numel(siteNames)
    name = siteNames{k};
    f = fullfile(csvDir, sprintf('dynamic_%s.mat', name));
    if ~isfile(f)
        warning('summarise_sites:missing', ...
            'No saved result for %s at %s -- skipped.', name, f);
        continue
    end
    S = load(f, 'OUTs');
    OUT = S.OUTs;
    A = OUT.apv.annual;
    C = OUT.compare.dynamic;

    rows(end+1,:) = { name, ...
        C.PV_total_MWh_ha, C.biogas_total_MWh_ha, ...
        C.LER_crop_PARstyle, C.LER_crop, C.LER_PV, C.LER_biogas, C.eLER, ...
        A.T_dig_mean_degC, 100*A.frac_time_mesophilic, ...
        100*A.heat_coverage, 100*A.pv_share_of_heat, ...
        100*A.frac_time_feed_limited, 100*A.frac_time_inhibited, ...
        A.pH_min, 100*A.x_CH4_mean, A.fruit_fresh_t_ha, ...
        100*A.periodicity.worst_rel, logical(A.periodicity.converged) }; %#ok<AGROW>
end

if isempty(rows)
    error('summarise_sites:noResults', ...
        'No saved results found in %s. Run run_all_sites first.', csvDir);
end

T = cell2table(rows, 'VariableNames', { ...
    'site','PV_MWh_ha','biogas_MWh_ha', ...
    'LER_crop_PAR','LER_crop_biomass','LER_PV','LER_biogas','eLER', ...
    'T_dig_mean_C','pct_time_mesophilic','pct_heat_covered', ...
    'pct_heat_from_PV','pct_time_feed_limited','pct_time_VFA_inhibited', ...
    'pH_min','pct_CH4','fruit_t_ha','pct_worst_drift','periodic'});

f = fullfile(csvDir, 'dynamic_all_sites.csv');
writetable(T, f);

fprintf('\n================ cross-site summary ================\n');
disp(T);
fprintf('saved: %s\n', f);

% CELL2TABLE turns a column of scalar logicals into a logical column, not a
% cell column, so index it directly -- brace indexing fails here.
per = T.periodic;
if iscell(per), per = [per{:}]; end
bad = T.site(~logical(per));
if ~isempty(bad)
    fprintf(['\nWARNING: not at periodic steady state: %s\n' ...
             'Increase opts.n_years before quoting their annual ' ...
             'figures.\n'], strjoin(string(bad), ', '));
end
end
