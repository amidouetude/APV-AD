function S = run_all_sites(siteNames, opts)
%RUN_ALL_SITES Run the dynamic model across the four study sites.
%
%   S = RUN_ALL_SITES() runs Konya, Almeria, Ouagadougou and Freiburg with
%   each site's DE optimum, then builds the cross-site comparison table.
%
%   S = RUN_ALL_SITES(SITENAMES, OPTS) restricts the sites or overrides the
%   options passed to MAIN_APVAD_DYNAMIC.
%
%   Each site's full result is saved by MAIN_APVAD_DYNAMIC before the next
%   one starts, so a failure late in the sequence never costs the earlier
%   sites' work; SUMMARISE_SITES can rebuild the table from what was saved.
%
%   A note on the temperature model. 'arrhenius' is the default everywhere
%   because it matches the static model's own temperature law and keeps the
%   benchmark clean. At Ouagadougou that choice deserves scrutiny: the
%   static model already places the digester at 33 degC on annual average,
%   so excursions past the mesophilic optimum are plausible, and the
%   Arrhenius form -- capped at one, with no descending limb -- cannot
%   penalise overheating. Re-run that site with
%       run_all_sites({'Ouagadougou'}, struct('temp_model','ctmi'))
%   and report both, rather than quietly picking one.
%
%   Runtime is roughly three to four minutes per site.
%
%   See also MAIN_APVAD_DYNAMIC, SUMMARISE_SITES, COMPARE_WITH_STATIC.

thisDir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(thisDir));

if nargin < 1 || isempty(siteNames)
    siteNames = {'Konya','Almeria','Ouagadougou','Freiburg'};
end
if nargin < 2, opts = struct(); end
if ~isfield(opts,'plots'),        opts.plots = false;       end
if ~isfield(opts,'save_results'), opts.save_results = true; end
if ~isfield(opts,'n_years'),      opts.n_years = 2;         end

S.results = struct();
ok = {};

for k = 1:numel(siteNames)
    name = siteNames{k};
    fprintf('\n########## %s (%d/%d) ##########\n', name, k, numel(siteNames));
    try
        OUT = main_apvad_dynamic(name, opts);
        S.results.(name) = OUT.compare;
        ok{end+1} = name; %#ok<AGROW>
    catch ME
        fprintf('FAILED %s: %s\n', name, ME.message);
        S.results.(name) = ME;
    end
end

if isempty(ok)
    S.table = table();
    return
end

S.table = summarise_sites(ok);
end
