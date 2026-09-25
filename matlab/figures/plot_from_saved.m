function OUT = plot_from_saved(siteName, csvDir, figDir)
%PLOT_FROM_SAVED Redraw the figures from a saved simulation result.
%
%   OUT = PLOT_FROM_SAVED(SITENAME) loads
%   outputs/csv/dynamic_<SITENAME>.mat, redraws both figures and returns
%   the stored result struct.
%
%   OUT = PLOT_FROM_SAVED(SITENAME, CSVDIR, FIGDIR) uses explicit folders.
%
%   A full-year run of the coupled model takes minutes; redrawing takes
%   less than a second. Separating the two means figure layout, colours and
%   detail windows can be iterated freely without ever re-integrating, and
%   it is also how figures are produced when the simulation itself ran on a
%   background worker with no display attached.
%
%   See also MAIN_APVAD_DYNAMIC, PLOT_APVAD_DYNAMICS.

thisDir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(thisDir));
projRoot = fileparts(thisDir);

if nargin < 2 || isempty(csvDir)
    csvDir = fullfile(projRoot, 'outputs', 'csv');
end
if nargin < 3 || isempty(figDir)
    figDir = fullfile(projRoot, 'outputs', 'figures');
end

f = fullfile(csvDir, sprintf('dynamic_%s.mat', siteName));
if ~isfile(f)
    error('plot_from_saved:notFound', ...
        'No saved result at %s. Run main_apvad_dynamic(''%s'') first.', ...
        f, siteName);
end

S = load(f, 'OUTs');
OUT = S.OUTs;
plot_apvad_dynamics(OUT, figDir);
end
