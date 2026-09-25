function poa = poa_irradiance(GHI, DNI, DHI, dec_deg, ha_deg, F_shad, ...
                              beta_deg, lat_deg, P)
%POA_IRRADIANCE Plane-of-array irradiance for a fixed due-south bifacial row.
%
%   POA = POA_IRRADIANCE(GHI, DNI, DHI, DEC_DEG, HA_DEG, F_SHAD, BETA_DEG,
%   LAT_DEG, P) returns a struct with fields G_front, G_rear, G_eff and
%   cos_theta_i (all element-wise on the input vectors).
%
%   Front surface, isotropic sky model:
%       G_front = DNI*cos(theta_i)
%               + DHI*(1 + cos(beta))/2
%               + GHI*rho_g*(1 - cos(beta))/2
%
%   Angle of incidence on a fixed, due-south collector (Duffie & Beckman):
%       cos(theta_i) = sin(delta)*sin(phi - beta)
%                    + cos(delta)*cos(H)*cos(phi - beta)
%
%   This exact closed form replaces the solar-noon-only approximation
%   sin(alpha_s + beta) used in the original Python code, which silently
%   assumed the sun stayed due south all day and over-predicted beam
%   irradiance every hour away from noon (see docs/poa_azimuth_fix_notes.md).
%
%   Rear surface, ground-reflected and reduced by the shaded ground
%   fraction, with the bifaciality factor applied in G_eff:
%       G_rear = GHI * rho_g * (1 - F_shad) * phi_rear
%       G_eff  = G_front + bifaciality * G_rear
%
%   See also SOLAR_ANGLES, SHADING_FRACTION, PV_MODULE.

beta_r = deg2rad(beta_deg);
lat_r  = deg2rad(lat_deg);
dec_r  = deg2rad(dec_deg);
ha_r   = deg2rad(ha_deg);

cos_theta_i = sin(dec_r).*sin(lat_r - beta_r) ...
            + cos(dec_r).*cos(ha_r).*cos(lat_r - beta_r);
cos_theta_i = min(max(cos_theta_i, 0), 1);

G_beam    = DNI .* cos_theta_i;
G_diffuse = DHI .* (1 + cos(beta_r))/2;
G_reflect = GHI .* P.rho_ground .* (1 - cos(beta_r))/2;
G_front   = max(G_beam + G_diffuse + G_reflect, 0);

G_rear = max(GHI .* P.rho_ground .* (1 - F_shad) .* P.phi_rear, 0);

poa.cos_theta_i = cos_theta_i;
poa.G_front = G_front;
poa.G_rear  = G_rear;
poa.G_eff   = G_front + P.bifaciality * G_rear;
end
