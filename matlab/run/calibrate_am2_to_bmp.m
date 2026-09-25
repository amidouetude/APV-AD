function [D, info] = calibrate_am2_to_bmp(D, P, HRT_days, c_VS_in_gL, BMP_target)
%CALIBRATE_AM2_TO_BMP Scale the methane yield coefficient k6 to a BMP assay.
%
%   [D, INFO] = CALIBRATE_AM2_TO_BMP(D, P, HRT_DAYS, C_VS_IN_GL) rescales
%   D.am2.k6 so that the AM2 digester, held at the reference temperature
%   with a steady feed, delivers the site's reference biochemical methane
%   potential P.BMP_ref (300 NmL CH4/g VS at 37 degC).
%
%   [D, INFO] = CALIBRATE_AM2_TO_BMP(..., BMP_TARGET) uses a different
%   target, for instance a measured BMP for tomato residue co-digested with
%   cattle manure.
%
%   Why calibrate here and not against the annual result
%   ---------------------------------------------------
%   The AM2 constants shipped in APVAD_DYN_PARAMS were identified on wine
%   distillery vinasse, so their absolute methane yield has no claim to be
%   right for tomato residue and manure. Something must set the scale. The
%   temptation is to tune the dynamic model until its annual methane output
%   matches the static model's -- which would be circular, and would throw
%   away the only interesting result by construction.
%
%   Instead the scale is set at the condition where BMP_ref is actually
%   DEFINED: a batch assay at 37 degC with substrate in excess of
%   inhibition. Everything downstream -- what the digester delivers at
%   Konya in February, how much of the year it spends below the mesophilic
%   band -- is then a prediction of the dynamic model, and its agreement or
%   disagreement with the static annual figure is a finding rather than a
%   fitted parameter.
%
%   k6 enters only the methane output equation qM = k6*mu2*X2 and none of
%   the mass balances, so the calibration is exactly linear and needs no
%   iteration:
%       k6_new = k6_old * BMP_target / BMP_model
%
%   INFO reports BMP before and after, the scale factor, and the steady
%   state used.
%
%   See also AM2_STEADY_STATE, APVAD_DYN_PARAMS.

if nargin < 5 || isempty(BMP_target)
    BMP_target = P.BMP_ref;
end

T_ref = D.temp.T_ref;
Dil   = 1 / HRT_days;

% Influent composition at the assay condition, built exactly as
% FEED_SCHEDULER builds it during a simulation.
S1in = D.feed.f_S1 * D.feed.COD_per_VS * c_VS_in_gL;
S2in = D.feed.S2in_mmolL;

ss = am2_steady_state(Dil, S1in, S2in, T_ref, D, P);
if ss.washout
    error('calibrate_am2_to_bmp:washout', ...
        ['The digester washes out at HRT = %.1f d with an influent of ' ...
         '%.2f g VS/L, so no steady methane yield exists to calibrate ' ...
         'against. Check the design variables.'], HRT_days, c_VS_in_gL);
end

% Specific methane yield at steady state.
%   qM [mmol/(L d)] * 22.414e-3 [Nm3/mol] * 1e3  -> Nm3 CH4 per m3 per day
%   fed VS         = Dil * c_VS_in                  kg VS per m3 per day
%   1 Nm3/kg VS    = 1000 NmL/g VS
molar_vol_Nm3 = 22.414e-3;
CH4_Nm3_m3_d  = ss.qM * molar_vol_Nm3;             % Nm3/(m3 d)
VS_fed_kg_m3_d = Dil * c_VS_in_gL;                 % kg VS/(m3 d)
BMP_model = 1000 * CH4_Nm3_m3_d / max(VS_fed_kg_m3_d, eps);   % NmL/g VS

scale = BMP_target / max(BMP_model, eps);

info.k6_before   = D.am2.k6;
info.BMP_model   = BMP_model;
info.BMP_target  = BMP_target;
info.scale       = scale;
info.steady_state= ss;
info.HRT_days    = HRT_days;
info.c_VS_in_gL  = c_VS_in_gL;

D.am2.k6 = D.am2.k6 * scale;
info.k6_after = D.am2.k6;

% Verify by recomputing (k6 does not feed back into the balances).
ss2 = am2_steady_state(Dil, S1in, S2in, T_ref, D, P);
info.BMP_after = 1000 * ss2.qM * molar_vol_Nm3 / max(VS_fed_kg_m3_d, eps);
end
