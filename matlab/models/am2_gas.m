function g = am2_gas(S2, Z, C, mu2, X2, V_dig_m3, D, P)
%AM2_GAS Gas-phase outputs of the AM2 digester: CH4, CO2, composition, pH.
%
%   G = AM2_GAS(S2, Z, C, MU2, X2, V_DIG_M3, D, P) returns a struct with
%
%     G.qM_mmol_L_d    specific methane production rate (mmol/(L d))
%     G.qC_mmol_L_d    specific carbon dioxide transfer rate (mmol/(L d))
%     G.Q_CH4_Nm3_d    methane volumetric flow (Nm3/d)
%     G.Q_CO2_Nm3_d    carbon dioxide volumetric flow (Nm3/d)
%     G.Q_biogas_Nm3_d total biogas flow (Nm3/d)
%     G.x_CH4          methane volume fraction (-)
%     G.P_CO2_atm      CO2 partial pressure in the headspace (atm)
%     G.CO2_mmol_L     dissolved inorganic carbon available as CO2
%     G.pH             digester pH (-)
%     G.E_CH4_kWh_d    methane energy production rate (kWh/d)
%
%   Methane leaves the liquid phase as fast as it is produced (its
%   solubility is negligible), so
%       qM = k6 * mu2 * X2.
%
%   Carbon dioxide is distributed between phases. AM2 closes the acid-base
%   system algebraically. Alkalinity Z is carried by bicarbonate and the
%   dissociated fatty acids, Z = [HCO3-] + S2, and total inorganic carbon
%   is C = [CO2] + [HCO3-]; eliminating bicarbonate gives
%       CO2   = C + S2 - Z
%       B     = Z - S2                     (bicarbonate)
%       phi   = CO2 + KH*Pt + qM/kLa
%       P_CO2 = ( phi - sqrt(phi^2 - 4*KH*Pt*CO2) ) / (2*KH)
%       qC    = kLa * (CO2 - KH*P_CO2)
%       pH    = -log10( Kb * CO2 / B )
%
%   Note the pH denominator is the BICARBONATE Z - S2, not C - Z. With a
%   well-buffered manure feed the two happen to be numerically close when
%   Cin ~ Zin, so the error is easy to miss -- until the influent
%   alkalinity is changed and pH runs to physically impossible values.
%
%   The discriminant is non-negative for all physical states: with qM = 0 it
%   reduces to (CO2 - KH*Pt)^2, and qM > 0 only increases it. It is clamped
%   at zero anyway so a rejected solver step cannot produce a complex value.
%
%   This is what makes the biogas methane fraction an OUTPUT of the dynamic
%   model rather than the fixed 0.60 that PARAMS.CH4_fraction imposes on the
%   static model, and it yields pH -- the variable an operator actually
%   watches -- for free.
%
%   Reference
%     Bernard, O. et al. (2001). Biotechnol. Bioeng., 75(4), 424-438.
%
%   See also AM2_KINETICS, APVAD_RHS.

a = D.am2;

% ---- Methane ---------------------------------------------------------
qM = a.k6 * mu2 .* max(X2, 0);                 % mmol/(L d)

% ---- Carbon dioxide / acid-base equilibrium --------------------------
CO2 = max(C + S2 - Z, a.c_floor);              % mmol/L

phi  = CO2 + a.KH*a.Pt + qM/a.kLa;
disc = max(phi.^2 - 4*a.KH*a.Pt*CO2, 0);
P_CO2 = (phi - sqrt(disc)) / (2*a.KH);         % atm
P_CO2 = max(P_CO2, 0);

qC = a.kLa * (CO2 - a.KH*P_CO2);               % mmol/(L d)
qC = max(qC, 0);

% ---- pH --------------------------------------------------------------
bicarb = max(Z - S2, a.c_floor);               % mmol/L
pH = -log10(max(a.Kb .* CO2 ./ bicarb, realmin));

% ---- Volumetric flows ------------------------------------------------
% mmol/(L d) * (V m3 * 1000 L/m3) / 1000 -> mol/d ; * 22.414e-3 -> Nm3/d
molar_vol_Nm3 = 22.414e-3;                     % Nm3/mol at 0 degC, 1 atm
Q_CH4 = qM * V_dig_m3 * molar_vol_Nm3;
Q_CO2 = qC * V_dig_m3 * molar_vol_Nm3;
Q_tot = Q_CH4 + Q_CO2;

x_CH4 = Q_CH4 ./ max(Q_tot, realmin);

% ---- Energy ----------------------------------------------------------
E_CH4_kWh_d = Q_CH4 * P.LHV_CH4_MJ_m3 / 3.6;   % MJ/d -> kWh/d

g.qM_mmol_L_d    = qM;
g.qC_mmol_L_d    = qC;
g.Q_CH4_Nm3_d    = Q_CH4;
g.Q_CO2_Nm3_d    = Q_CO2;
g.Q_biogas_Nm3_d = Q_tot;
g.x_CH4          = x_CH4;
g.P_CO2_atm      = P_CO2;
g.CO2_mmol_L     = CO2;
g.pH             = pH;
g.E_CH4_kWh_d    = E_CH4_kWh_d;
end
