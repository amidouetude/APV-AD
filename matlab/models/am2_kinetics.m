function k = am2_kinetics(S1, S2, T_dig, D, P)
%AM2_KINETICS Specific growth rates of the two AM2 bacterial populations.
%
%   K = AM2_KINETICS(S1, S2, T_DIG, D, P) returns a struct with fields
%   mu1, mu2 (1/d) and fT (-), the temperature factor applied to both.
%
%   Acidogenesis -- Monod, no inhibition:
%       mu1 = fT * mu1max * S1 / (KS1 + S1)
%
%   Methanogenesis -- Haldane, substrate-inhibited by VFA:
%       mu2 = fT * mu2max * S2 / (KS2 + S2 + S2^2/KI2)
%
%   The Haldane term is the reason a dynamic model earns its keep here.
%   The static annual balance cannot represent VFA accumulation at all: it
%   converts volatile solids to methane through a single temperature-
%   corrected BMP, so an overload simply produces proportionally more gas.
%   In the dynamic model an organic overload raises S2, mu2 falls once S2
%   passes sqrt(KS2*KI2), less VFA is consumed, S2 rises further, and the
%   digester can fall into the acidified state that operators actually
%   observe. Souring is a trajectory, not an annual average.
%
%   Concentrations are floored at D.am2.c_floor before use so that the
%   right-hand side remains smooth and finite if the solver probes slightly
%   negative values during a rejected step.
%
%   Reference
%     Bernard, O. et al. (2001). Biotechnol. Bioeng., 75(4), 424-438.
%
%   See also GROWTH_TEMP_FACTOR, AM2_GAS, APVAD_RHS.

a  = D.am2;
S1 = max(S1, a.c_floor);
S2 = max(S2, a.c_floor);

fT = growth_temp_factor(T_dig, D, P);

mu1 = fT .* a.mu1max .* S1 ./ (a.KS1 + S1);
mu2 = fT .* a.mu2max .* S2 ./ (a.KS2 + S2 + (S2.^2)/a.KI2);

k.mu1 = mu1;
k.mu2 = mu2;
k.fT  = fT;
end
