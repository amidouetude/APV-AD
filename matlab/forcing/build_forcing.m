function F = build_forcing(clim, site, design, P, D, mode)
%BUILD_FORCING Precompute exogenous signals and wrap them as interpolants.
%
%   F = BUILD_FORCING(CLIM, SITE, DESIGN, P, D)          builds the APV case
%   F = BUILD_FORCING(CLIM, SITE, DESIGN, P, D, 'reference')  builds the
%       open-field reference case (no inter-row shading, dense-packed PV),
%       which supplies the denominators of every LER component.
%
%   Inputs
%     CLIM    struct from LOAD_CLIMATE
%     SITE    struct from APVAD_SITES(name)
%     DESIGN  struct with fields beta_deg, d_row_m, H_m_m, V_dig_m3,
%             HRT_days, f_PV_heat
%     P, D    APVAD_PARAMS and APVAD_DYN_PARAMS
%
%   Everything that depends on time but not on the model states is
%   evaluated once here, on the hourly climate grid, and exposed to the
%   solver as griddedInterpolant objects. The stiff solver may call the
%   right-hand side tens of thousands of times per simulated year, so doing
%   the solar geometry, transposition and ET0 algebra inside the RHS would
%   dominate the runtime for no benefit -- none of it depends on the state.
%
%   Interpolants are linear in time with nearest-neighbour extrapolation,
%   so the RHS stays defined if the solver probes slightly outside the
%   climate record.
%
%   Output fields (all interpolants are functions of t in days, t = 0 at
%   1 January 00:00)
%     F.T_amb, F.WS, F.GHI, F.G_eff, F.F_shad, F.PAR_AV, F.PAR_open
%     F.ET0        mm/d   reference ET rate
%     F.rain       mm/d   precipitation rate
%     F.Kc         -      FAO-56 crop coefficient
%     F.in_season  0..1   growing-season indicator
%     F.T_day      degC   daily mean air temperature
%     F.harvest    1/d    unit-area residue transfer pulse
%     F.n_modules, F.GCR, F.A_pv_m2, F.mode, F.doy_plant, F.doy_harvest
%     F.manure_VS_kg_d, F.c_VS_in_gL, F.Q_in_m3_d, F.A_wall_m2, F.V_dig_m3
%
%   See also LOAD_CLIMATE, SOLAR_ANGLES, POA_IRRADIANCE, ET0_HARGREAVES.

if nargin < 6 || isempty(mode)
    mode = 'apv';
end
isRef = strcmpi(mode, 'reference');

t   = clim.t_day(:);
n   = numel(t);
hh  = hour(clim.time) + minute(clim.time)/60;
doy = double(clim.doy(:));

% ---------------------------------------------------------------------
% Solar geometry and shading
% ---------------------------------------------------------------------
sun = solar_angles(doy, hh, site.lat, site.lon);

beta = design.beta_deg;
if isRef
    % Dense-packed open-field reference: rows spaced at exactly the module
    % projected width, so modules never shade the ground between rows.
    d_row  = max(P.module_length_m*cos(deg2rad(beta)), 0.5);
    F_shad = zeros(n,1);
    GCR    = P.module_length_m*cos(deg2rad(beta))/d_row;
else
    d_row = design.d_row_m;
    [F_shad, GCR] = shading_fraction(sun.elevation_deg, beta, d_row, ...
                                     design.H_m_m, P);
end

n_modules = modules_per_ha(d_row, P);

% ---------------------------------------------------------------------
% Plane-of-array irradiance
% ---------------------------------------------------------------------
poa = poa_irradiance(clim.GHI, clim.DNI, clim.DHI, ...
                     sun.declination_deg, sun.hour_angle_deg, F_shad, ...
                     beta, site.lat, P);

% ---------------------------------------------------------------------
% PAR and reference evapotranspiration
% ---------------------------------------------------------------------
PAR_open = clim.GHI * P.f_PAR;              % W/m2
PAR_AV   = PAR_open .* (1 - F_shad);        % W/m2 beneath the array

[ET0_h, ~, ~] = et0_hargreaves(clim.time, clim.GHI, clim.T2m, site.lat);
ET0_rate = ET0_h * 24;                      % mm/h -> mm/d
rain_rate = clim.P_mm * 24;                 % mm/h -> mm/d

% ---------------------------------------------------------------------
% Daily mean air temperature (drives thermal time)
% ---------------------------------------------------------------------
[~, ~, dayIdx] = unique(dateshift(clim.time, 'start', 'day'));
Tday = accumarray(dayIdx, clim.T2m, [], @mean);
T_day_h = Tday(dayIdx);

% ---------------------------------------------------------------------
% Season calendar, FAO-56 crop coefficient, harvest pulse
% ---------------------------------------------------------------------
yr = year(clim.time(1));
doy_plant   = day(datetime(yr, str2double(site.planting_date(1:2)), ...
                                str2double(site.planting_date(4:5))), 'dayofyear');
doy_harvest = day(datetime(yr, str2double(site.harvest_date(1:2)), ...
                                str2double(site.harvest_date(4:5))), 'dayofyear');
seasonLen = doy_harvest - doy_plant;
if seasonLen <= 0
    seasonLen = seasonLen + 365;            % season wraps the new year
end

frac = (doy - doy_plant) / seasonLen;       % 0 at planting, 1 at harvest
inSeason = double(clim.in_season(:));

% FAO-56 four-stage single crop coefficient, linear on the season fraction
Kc = zeros(n,1);
Kc(frac >= 0    & frac < 0.20) = D.water.Kc_ini;
mid = frac >= 0.20 & frac < 0.45;
Kc(mid) = D.water.Kc_ini + (D.water.Kc_mid - D.water.Kc_ini) ...
          .* (frac(mid) - 0.20)/0.25;
Kc(frac >= 0.45 & frac < 0.75) = D.water.Kc_mid;
lat_ = frac >= 0.75 & frac <= 1.0;
Kc(lat_) = D.water.Kc_mid + (D.water.Kc_end - D.water.Kc_mid) ...
           .* (frac(lat_) - 0.75)/0.25;
Kc = Kc .* inSeason;

% Harvest removal-rate constant: a Gaussian centred on the harvest date,
% scaled so its time-integral equals D.feed.harvest_strength. Removal is
% first-order in standing biomass, so that integral sets the fraction
% harvested, 1 - exp(-strength).
%
% Spreading the harvest over a window rather than applying a discontinuous
% jump is not cosmetic: a variable-step stiff solver would either step
% straight over an instantaneous transfer or collapse its step size trying
% to resolve it.
sw   = D.feed.harvest_window_d/2;
gpul = exp(-0.5*((doy - doy_harvest)/sw).^2);
area = trapz(t, gpul);
if area <= 0
    harvest_pulse = zeros(n,1);
else
    harvest_pulse = D.feed.harvest_strength * gpul/area;   % 1/d
end

% ---------------------------------------------------------------------
% Digester sizing quantities
% ---------------------------------------------------------------------
V_dig = design.V_dig_m3;
Q_in  = V_dig / design.HRT_days;            % m3/d

% Cylindrical digester with height = diameter, as in the static model
D_m    = (4*V_dig/pi)^(1/3);
A_wall = pi*D_m*D_m + 2*pi*(D_m/2)^2;       % lateral + top + bottom, m2

% Annual VS availability, used to set the feed dilution so that the
% nominal organic loading rate matches the static model's
% OLR = VS_annual / (365 * V_dig).
%
% The default basis is the site's reported open-field residue yield. Under
% an array the crop yields less, so the design basis and the realised feed
% disagree, and the store simply runs dry for part of the year. MAIN_APVAD_
% DYNAMIC therefore runs the agrivoltaic case twice: a first pass on the
% nominal basis to discover the realised residue production, then the
% reported run with VS_ann_kg set to that realised figure, so that the
% dilution, the loading rate and the feedstock actually available are
% mutually consistent -- exactly the design a farmer would end up with.
VS_res_nom = site.fruit_yield_t_ha*1000 * site.R_res * site.DM * site.VS_fraction;
VS_man_nom = site.manure_input_t_ha*1000 * site.manure_DM * site.manure_VS_fraction;

if isfield(design, 'VS_ann_kg') && ~isempty(design.VS_ann_kg)
    VS_ann_nom = design.VS_ann_kg;
else
    VS_ann_nom = VS_res_nom + VS_man_nom;                 % kg VS/ha/yr
end
c_VS_in    = VS_ann_nom / (365 * max(Q_in, eps));         % kg/m3 = g/L

manure_VS_kg_d = VS_man_nom / 365;                        % kg VS/d

% ---------------------------------------------------------------------
% Heater sizing (see APVAD_DYN_PARAMS section 4)
% ---------------------------------------------------------------------
% Design load = wall loss at the coldest monthly mean ambient temperature
% + sensible heat to raise the influent from T_in to the setpoint.
[~, ~, monIdx] = unique(month(clim.time));
T_month_mean   = accumarray(monIdx, clim.T2m, [], @mean);
T_amb_design   = min(T_month_mean);

U_eff = P.U_eff_W_m2K * D.thermal.insulation_scale;
Q_loss_design = U_eff * A_wall * (D.thermal.T_set_degC - T_amb_design);   % W
Q_flow_design = Q_in * P.rho_substrate * P.Cp_substrate ...
                * max(D.thermal.T_set_degC - D.feed.T_in_degC, 0) / 86.4; % W
Q_design = max(Q_loss_design + Q_flow_design, 1);

if isempty(D.thermal.Q_max_W)
    Q_max_W = D.thermal.Q_max_factor * Q_design;
else
    Q_max_W = D.thermal.Q_max_W;
end
if isempty(D.thermal.Kp_W_per_K)
    Kp_W_per_K = Q_max_W / max(D.thermal.prop_band_K, eps);
else
    Kp_W_per_K = D.thermal.Kp_W_per_K;
end

% ---------------------------------------------------------------------
% Wrap as interpolants
% ---------------------------------------------------------------------
mk = @(v) griddedInterpolant(t, double(v(:)), 'linear', 'nearest');

F.t_day      = t;
F.T_amb      = mk(clim.T2m);
F.WS         = mk(clim.WS);
F.GHI        = mk(clim.GHI);
F.G_eff      = mk(poa.G_eff);
F.G_front    = mk(poa.G_front);
F.F_shad     = mk(F_shad);
F.PAR_AV     = mk(PAR_AV);
F.PAR_open   = mk(PAR_open);
F.ET0        = mk(ET0_rate);
F.rain       = mk(rain_rate);
F.Kc         = mk(Kc);
F.in_season  = mk(inSeason);
F.T_day      = mk(T_day_h);
F.harvest    = mk(harvest_pulse);

F.mode         = mode;
F.site         = site;
F.design       = design;
F.n_modules    = n_modules;
F.GCR          = GCR;
F.d_row_m      = d_row;
F.A_pv_m2      = n_modules * P.module_length_m * P.module_width_m;
F.doy_plant    = doy_plant;
F.doy_harvest  = doy_harvest;
F.season_len_d = seasonLen;
F.V_dig_m3     = V_dig;
F.Q_in_m3_d    = Q_in;
F.A_wall_m2    = A_wall;
F.c_VS_in_gL   = c_VS_in;
F.manure_VS_kg_d = manure_VS_kg_d;
F.VS_res_nom_kg  = VS_res_nom;
F.Q_design_W     = Q_design;
F.Q_max_W        = Q_max_W;
F.Kp_W_per_K     = Kp_W_per_K;
F.T_amb_design_degC = T_amb_design;
F.VS_ann_nom_kg  = VS_ann_nom;
F.t_end_day      = t(end);

% Period used to wrap simulated time back onto the climate record, so that
% a multi-year run repeats the same TMY year. A typical meteorological year
% is by construction a representative year, not a trajectory, so repeating
% it is the right way to reach a periodic steady state -- the alternative,
% integrating one year from arbitrary initial conditions, reports the
% digester's start-up transient as if it were an annual result.
F.t_period = t(end) + median(diff(t));
end

% =====================================================================
function n = modules_per_ha(d_row_m, P)
%MODULES_PER_HA Module count on one hectare, identical to the Python model.
%   Rows run east-west; n_rows = field_depth/d_row, n_cols = row_length /
%   module_length.
field_depth = P.land_area_m2 / P.row_length_m;
n_rows = floor(field_depth / d_row_m);
n_cols = floor(P.row_length_m / P.module_length_m);
n = max(1, n_rows*n_cols);
end
