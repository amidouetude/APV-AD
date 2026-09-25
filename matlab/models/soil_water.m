function w = soil_water(Dr, ET0, Kc, F_shad, rain, inSeason, D, P)
%SOIL_WATER FAO-56 single-coefficient root-zone water balance.
%
%   W = SOIL_WATER(DR, ET0, KC, F_SHAD, RAIN, INSEASON, D, P) returns
%
%     W.dDr_dt      root-zone depletion derivative (mm/d)
%     W.Ks          water-stress coefficient (-)
%     W.ETc_mm_d    actual crop evapotranspiration under the array (mm/d)
%     W.ETc_open    the same crop without array shading (mm/d)
%     W.irr_mm_d    applied irrigation (mm/d)
%     W.perc_mm_d   deep percolation (mm/d)
%     W.TAW, W.RAW  total and readily available water (mm)
%
%   Depletion form of the balance (FAO-56 Ch. 8), where Dr is the water
%   deficit of the root zone relative to field capacity:
%
%       dDr/dt = ETc_adj - (rain + irrigation) + percolation
%       TAW    = 1000*(theta_fc - theta_wp)*Zr
%       RAW    = p*TAW
%       Ks     = 1                            for Dr <= RAW
%              = (TAW - Dr)/(TAW - RAW)       for Dr >  RAW
%       ETc_adj= Ks * Kc * ET0 * (1 - alpha_shade*F_shad)
%
%   The shading term is the same one the static model applies, but here it
%   acts on a state. That distinction is the point: the static model
%   reports water SAVED, an accounting quantity, because it never tracks
%   how much water is in the soil. Carrying Dr as a state means the shading
%   benefit shows up where it matters agronomically -- as fewer hours spent
%   below the stress threshold, and as irrigation events deferred or
%   skipped -- and it feeds back into growth through Ks.
%
%   Irrigation refills towards field capacity whenever depletion passes
%   RAW, at a bounded application rate, and only within the growing season.
%   Percolation removes any water above field capacity (Dr < 0) with a
%   short relaxation time, which keeps Dr >= 0 without a hard clamp.
%
%   Reference
%     Allen, R.G. et al. (1998). FAO Irrigation and Drainage Paper 56.
%
%   See also CROP_GROWTH, SMOOTH_STEP.

TAW = 1000*(D.water.theta_fc - D.water.theta_wp)*D.water.Zr_m;
RAW = D.water.p_depletion * TAW;

% ---- Water stress coefficient ----------------------------------------
Ks_lin = (TAW - Dr) / max(TAW - RAW, eps);
Ks_lin = min(max(Ks_lin, 0), 1);
% Blend between the unstressed value and the linear decline so that Ks is
% continuously differentiable at Dr = RAW.
sStress = smooth_step(Dr - RAW, max(0.02*TAW, eps));
Ks = (1 - sStress)*1 + sStress*Ks_lin;

% ---- Crop evapotranspiration -----------------------------------------
ETc_open = Ks * Kc * ET0;
ETc      = ETc_open * (1 - P.alpha_shade*F_shad);

% ---- Irrigation -------------------------------------------------------
if D.water.irrigate
    demand = smooth_step(Dr - RAW, max(0.05*TAW, eps));
    irr    = D.water.irr_rate_mm_d * demand * inSeason;
    % Never apply more than is needed to reach field capacity.
    irr    = irr .* smooth_step(Dr, max(0.02*TAW, eps));
else
    irr = 0;
end

% ---- Deep percolation -------------------------------------------------
tau_perc = 0.05;                          % d
perc = max(-Dr, 0) / tau_perc;

% ---- Balance ----------------------------------------------------------
dDr_dt = ETc - (rain + irr) + perc;

w.dDr_dt    = dDr_dt;
w.Ks        = Ks;
w.ETc_mm_d  = ETc;
w.ETc_open  = ETc_open;
w.irr_mm_d  = irr;
w.perc_mm_d = perc;
w.TAW       = TAW;
w.RAW       = RAW;
end
