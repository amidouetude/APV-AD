function [ET0_h, ET0_d, Ra_d] = et0_hargreaves(time, GHI, T2m, lat_deg)
%ET0_HARGREAVES Reference evapotranspiration, Hargreaves-Samani / FAO-56.
%
%   [ET0_H, ET0_D, RA_D] = ET0_HARGREAVES(TIME, GHI, T2M, LAT_DEG)
%
%       ET0_H  n x 1  hourly reference ET (mm/h)
%       ET0_D  n x 1  the day's total ET0 repeated on every hour (mm/d)
%       RA_D   n x 1  extraterrestrial radiation (MJ/m2/d), repeated hourly
%
%   Daily reference ET (FAO-56 Eq. 52):
%       ET0 = 0.0023 * 0.408 * Ra * (Tmean + 17.8) * sqrt(Tmax - Tmin)
%
%   Extraterrestrial radiation Ra is computed astronomically from latitude
%   and day of year (FAO-56 Eqs. 21-25), not inferred from measured GHI.
%   Both points matter and both were corrections made to the Python model:
%
%     * Ra was previously estimated as GHI/0.75, assuming a constant
%       clearness index. Actual clearness indices at these four sites are
%       0.52-0.63, so Ra -- and therefore ET0 -- was underestimated by
%       15-34 %, worst at the cloudiest site.
%     * The 0.408 factor converting MJ/m2/d to mm/d (the reciprocal of the
%       latent heat of vaporisation, 2.45 MJ/kg) was missing entirely.
%
%   See docs/et0_ra_fix_notes.md for the derivation and quantified impact.
%
%   The daily total is distributed over the day in proportion to GHI, so
%   ET occurs during daylight hours only.
%
%   Reference
%     Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998). Crop
%       evapotranspiration -- Guidelines for computing crop water
%       requirements. FAO Irrigation and Drainage Paper 56.

doy = day(time, 'dayofyear');
[~, ~, dayIdx] = unique(dateshift(time, 'start', 'day'));

% ---- Extraterrestrial radiation Ra (MJ/m2/d), FAO-56 Eqs. 21-25 -------
lat_r   = deg2rad(lat_deg);
dr      = 1 + 0.033*cos(2*pi*doy/365);
delta_r = 0.409*sin(2*pi*doy/365 - 1.39);
omega_s = acos(min(max(-tan(lat_r).*tan(delta_r), -1), 1));
G_sc    = 0.0820;                                   % MJ/(m2 min)

Ra_d = (24*60/pi) * G_sc * dr .* ( ...
        omega_s.*sin(lat_r).*sin(delta_r) ...
      + cos(lat_r).*cos(delta_r).*sin(omega_s));

% ---- Daily air-temperature statistics ---------------------------------
Tmean_d = accumarray(dayIdx, T2m, [], @mean);
Tmax_d  = accumarray(dayIdx, T2m, [], @max);
Tmin_d  = accumarray(dayIdx, T2m, [], @min);
Trange  = max(Tmax_d - Tmin_d, 0);

Tmean_h  = Tmean_d(dayIdx);
Trange_h = Trange(dayIdx);

% ---- Hargreaves-Samani daily ET0 (mm/d) -------------------------------
ET0_d = max(0.0023 * 0.408 * Ra_d .* (Tmean_h + 17.8) .* sqrt(Trange_h), 0);

% ---- Distribute over daylight hours in proportion to GHI --------------
GHI_d   = accumarray(dayIdx, max(GHI, 0), [], @sum);
GHI_d_h = max(GHI_d(dayIdx), 1);

ET0_h = zeros(size(GHI));
lit   = GHI > 1;
ET0_h(lit) = ET0_d(lit) .* GHI(lit) ./ GHI_d_h(lit);
end
