function clim = load_climate(siteName, csvDir)
%LOAD_CLIMATE Read one site's hourly climate file into a clean struct.
%
%   CLIM = LOAD_CLIMATE(SITENAME) reads outputs/csv/hourly_<SITENAME>.csv,
%   the file produced by the existing 01_load_data.ipynb notebook, from the
%   project root inferred from this file's location.
%
%   CLIM = LOAD_CLIMATE(SITENAME, CSVDIR) reads from an explicit directory.
%
%   The returned struct has fields
%       t_day       n x 1  time in days since 1 January 00:00 (t = 0)
%       time        n x 1  datetime
%       T2m         n x 1  air temperature (degC)
%       GHI,DNI,DHI n x 1  irradiance components (W/m2)
%       WS          n x 1  wind speed (m/s)
%       RH          n x 1  relative humidity (%) if present
%       P_mm        n x 1  precipitation (mm/h), zeros if not present
%       in_season   n x 1  logical growing-season mask
%       doy         n x 1  day of year
%       dt_h        scalar sampling interval in hours
%       site        char   site name
%
%   PVGIS TMY files carry no precipitation column. When one is absent
%   P_mm is filled with zeros and a warning is issued: the soil-water
%   balance is then irrigation-driven only, which must be stated as a
%   limitation for humid sites (Freiburg in particular).
%
%   The typical meteorological year is spliced month by month from
%   DIFFERENT calendar years -- hourly_Konya.csv draws on 2005, 2008, 2009,
%   2010, 2012, 2013, 2017 and 2018 -- so its raw timestamps jump backwards
%   at every month boundary and are useless as a time axis. Worse, months
%   taken from a leap year carry a day-of-year offset by one, which is why
%   the file's own day_of_year column contains 366 values but only 364
%   distinct ones, with 24 duplicated hours.
%
%   This loader therefore rebuilds the calendar by mapping each row's
%   (month, day, hour) onto a single non-leap reference year and sorting.
%   The result is a strictly increasing 365-day axis, which is what the
%   interpolants and the ODE solver require. Any 29 February row is
%   dropped, since a typical year has no such date.
%
%   See also BUILD_FORCING.

arguments
    siteName (1,:) char
    csvDir   (1,:) char = ''
end

if isempty(csvDir)
    here    = fileparts(fileparts(mfilename('fullpath')));  % matlab/
    root    = fileparts(here);                              % project root
    csvDir  = fullfile(root, 'outputs', 'csv');
end

f = fullfile(csvDir, sprintf('hourly_%s.csv', siteName));
if ~isfile(f)
    error('load_climate:fileNotFound', ...
        ['Climate file not found:\n  %s\n' ...
         'Run 01_load_data.ipynb first, or pass an explicit csvDir.'], f);
end

opts = detectImportOptions(f, 'TextType', 'string');
opts = setvartype(opts, 'time_UTC', 'datetime');
if any(strcmp(opts.VariableNames, 'in_season'))
    opts = setvartype(opts, 'in_season', 'string');
end
T = readtable(f, opts);

% ---- Rebuild a monotone calendar on a non-leap reference year --------
REF_YEAR = 2001;                            % arbitrary non-leap year
tv = T.time_UTC;
isFeb29 = (month(tv) == 2) & (day(tv) == 29);
if any(isFeb29)
    T  = T(~isFeb29, :);
    tv = tv(~isFeb29);
end

t_ref = datetime(REF_YEAR, month(tv), day(tv)) + hours(hour(tv)) ...
        + minutes(minute(tv));
[t_ref, ord] = sort(t_ref);
T = T(ord, :);

if any(diff(t_ref) <= 0)
    error('load_climate:duplicateTimestamps', ...
        ['The reconstructed calendar for %s still contains duplicate or ' ...
         'out-of-order hours. Inspect the source file: a typical ' ...
         'meteorological year should hold exactly one row per hour.'], ...
        siteName);
end

clim.site  = siteName;
clim.time  = t_ref;
clim.t_day = days(t_ref - datetime(REF_YEAR,1,1));
clim.ref_year = REF_YEAR;
clim.source_years = unique(year(tv))';

req = {'T2m','GHI','DNI','DHI','WS'};
for k = 1:numel(req)
    if ~ismember(req{k}, T.Properties.VariableNames)
        error('load_climate:missingColumn', ...
              'Required column "%s" not found in %s.', req{k}, f);
    end
    clim.(req{k}) = double(T.(req{k}));
end

% DNI is stored as -0.0 at night in the PVGIS export; clip the negatives.
clim.GHI = max(clim.GHI, 0);
clim.DNI = max(clim.DNI, 0);
clim.DHI = max(clim.DHI, 0);
clim.WS  = max(clim.WS,  0);

if ismember('RH', T.Properties.VariableNames)
    clim.RH = double(T.RH);
else
    clim.RH = nan(size(clim.T2m));
end

% Precipitation -- optional.
rainCandidates = {'P_mm','precip_mm','precipitation','RR','P'};
clim.P_mm = zeros(size(clim.T2m));
found = false;
for k = 1:numel(rainCandidates)
    if ismember(rainCandidates{k}, T.Properties.VariableNames)
        clim.P_mm = double(T.(rainCandidates{k}));
        found = true;
        break
    end
end
if ~found
    warning('load_climate:noPrecipitation', ...
        ['No precipitation column in %s. The soil-water balance will be ' ...
         'irrigation-driven only; report this as a model limitation.'], ...
        sprintf('hourly_%s.csv', siteName));
end

% Growing-season mask, stored by the Python loader as the strings
% "True"/"False".
if ismember('in_season', T.Properties.VariableNames)
    v = T.in_season;
    if isstring(v) || iscellstr(v)
        clim.in_season = strcmpi(string(v), "True");
    elseif islogical(v)
        clim.in_season = v;
    else
        clim.in_season = double(v) > 0.5;
    end
else
    clim.in_season = false(size(clim.T2m));
end

clim.doy  = day(clim.time, 'dayofyear');
clim.dt_h = hours(median(diff(clim.time)));
end
