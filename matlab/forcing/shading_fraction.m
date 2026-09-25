function [F_shad, GCR] = shading_fraction(elev_deg, beta_deg, d_row_m, H_m_m, P)
%SHADING_FRACTION Ground shading fraction beneath an APV row array.
%
%   [F_SHAD, GCR] = SHADING_FRACTION(ELEV_DEG, BETA_DEG, D_ROW_M, H_M_M, P)
%
%       F_shad = min(1, H_m * cos(dGamma) / (d_row * tan(alpha_s)))
%       GCR    = module_length * cos(beta) / d_row
%
%   with dGamma = 0 for due-south rows. F_shad is zero at night
%   (alpha_s <= 0) and small at summer midday when the sun is high.
%
%   P is the struct returned by APVAD_PARAMS (module geometry).
%
%   This is the same expression as compute_shading() in the Python model,
%   evaluated here on the hourly grid before the ODE integration so that
%   the solver only ever interpolates a precomputed signal.
%
%   See also SOLAR_ANGLES, POA_IRRADIANCE, BUILD_FORCING.

delta_gamma = 0;                       % due-south rows
tan_alpha   = tan(deg2rad(min(max(elev_deg, 0.1), 90)));

F_raw  = H_m_m * cos(deg2rad(delta_gamma)) ./ (d_row_m * tan_alpha);
F_shad = min(max(F_raw, 0), 1);
F_shad(elev_deg <= 0) = 0;

GCR = (P.module_length_m * cos(deg2rad(beta_deg))) / d_row_m;
end
