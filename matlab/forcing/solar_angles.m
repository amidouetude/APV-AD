function sun = solar_angles(t_day, hour_utc, lat_deg, lon_deg)
%SOLAR_ANGLES Declination, hour angle and solar elevation.
%
%   SUN = SOLAR_ANGLES(T_DAY, HOUR_UTC, LAT_DEG, LON_DEG) returns a struct
%   with fields declination_deg, hour_angle_deg, elevation_deg and
%   eot_min, evaluated element-wise on the input time vectors.
%
%   T_DAY    day of year (1 = 1 January), any real vector
%   HOUR_UTC clock hour UTC, including fractional minutes
%
%   Declination follows Spencer (1971). The hour angle is computed from
%   TRUE LOCAL SOLAR TIME, i.e. clock time corrected for site longitude
%   (4 min per degree east) and for the equation of time. This is the
%   corrected formulation adopted in the Python model after the pvlib
%   cross-check reported in docs/pvlib_validation.md; using the raw UTC
%   hour instead implicitly places solar noon at 12:00 UTC everywhere,
%   which is only true on the Greenwich meridian and biased annual
%   front-surface POA by about -9.6 % at Konya.
%
%   References
%     Spencer, J.W. (1971). Fourier series representation of the position
%       of the sun. Search, 2(5), 172.
%     Duffie, J.A., Beckman, W.A. (2013). Solar Engineering of Thermal
%       Processes, 4th ed., Wiley.
%
%   See also POA_IRRADIANCE, SHADING_FRACTION.

B = 2*pi*(t_day - 1)/365;

dec = (180/pi) * ( ...
      0.006918 ...
    - 0.399912*cos(B)   + 0.070257*sin(B) ...
    - 0.006758*cos(2*B) + 0.000907*sin(2*B) ...
    - 0.002697*cos(3*B) + 0.001480*sin(3*B));

eot_min = 229.18 * ( ...
      0.000075 ...
    + 0.001868*cos(B)   - 0.032077*sin(B) ...
    - 0.014615*cos(2*B) - 0.040849*sin(2*B));

solar_time = hour_utc + lon_deg/15 + eot_min/60;   % true local solar time (h)
hour_angle = (solar_time - 12) * 15;               % deg, 0 at solar noon

lat_r = deg2rad(lat_deg);
dec_r = deg2rad(dec);
ha_r  = deg2rad(hour_angle);

sin_elev = sin(lat_r).*sin(dec_r) + cos(lat_r).*cos(dec_r).*cos(ha_r);
elev     = rad2deg(asin(min(max(sin_elev, -1), 1)));

sun.declination_deg = dec;
sun.eot_min         = eot_min;
sun.hour_angle_deg  = hour_angle;
sun.elevation_deg   = elev;
end
