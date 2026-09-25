function ss = am2_steady_state(Dil, S1in, S2in, T_degC, D, P)
%AM2_STEADY_STATE Closed-form steady state of the AM2 digester.
%
%   SS = AM2_STEADY_STATE(DIL, S1IN, S2IN, T_DEGC, D, P) returns the
%   non-washout equilibrium of the AM2 mass balances at constant dilution
%   rate DIL (1/d), constant influent and constant temperature, with fields
%   S1, X1, S2, X2, mu1, mu2, qM (mmol/(L d)) and a logical WASHOUT flag.
%
%   At equilibrium each population grows exactly as fast as it is removed:
%       mu1 = alpha*Dil + kd1,     mu2 = alpha*Dil + kd2
%
%   Inverting Monod for the acidogens is immediate:
%       S1* = KS1*mu1* / (fT*mu1max - mu1*)
%
%   Inverting Haldane for the methanogens gives a quadratic,
%       (m/KI2)*S2^2 + (m-1)*S2 + m*KS2 = 0,   m = mu2*/(fT*mu2max)
%   whose two roots are the two branches of the inhibition curve. The
%   SMALLER root is the healthy, locally stable operating point; the larger
%   root lies on the descending limb where more VFA means slower growth,
%   and is unstable. That pair of roots is the reactor's bistability, and
%   it is the structural reason a digester can sour irreversibly after a
%   transient overload -- something no annual-average model can exhibit.
%
%   Washout (no real positive root, or mu*max below the removal rate) is
%   reported rather than silently returning a complex or negative state.
%
%   This closed form is used to calibrate k6 against the laboratory BMP
%   assay, and by TESTAM2STEADYSTATEMATCHESINTEGRATION as an independent
%   check on the integrated model.
%
%   See also AM2_KINETICS, CALIBRATE_AM2_TO_BMP.

a  = D.am2;
fT = growth_temp_factor(T_degC, D, P);

mu1s = a.alpha*Dil + a.kd1;
mu2s = a.alpha*Dil + a.kd2;

ss.washout = false;

% ---- Acidogens: Monod inversion --------------------------------------
den1 = fT*a.mu1max - mu1s;
if den1 <= 0
    ss.washout = true;
    ss = fill_washout(ss, S1in, S2in);
    return
end
S1 = a.KS1 * mu1s / den1;
if S1 >= S1in
    ss.washout = true;
    ss = fill_washout(ss, S1in, S2in);
    return
end
X1 = Dil*(S1in - S1) / (a.k1*mu1s);

% ---- Methanogens: Haldane inversion ----------------------------------
m = mu2s / max(fT*a.mu2max, eps);
A = m/a.KI2;
B = m - 1;
C = m*a.KS2;

disc = B^2 - 4*A*C;
if disc < 0 || A <= 0
    ss.washout = true;
    ss = fill_washout(ss, S1in, S2in);
    return
end
r1 = (-B - sqrt(disc)) / (2*A);
r2 = (-B + sqrt(disc)) / (2*A);
roots_pos = sort([r1 r2]);
S2 = roots_pos(roots_pos > 0);
if isempty(S2)
    ss.washout = true;
    ss = fill_washout(ss, S1in, S2in);
    return
end
S2 = S2(1);                       % stable (lower) branch

X2 = (Dil*(S2in - S2) + a.k2*mu1s*X1) / (a.k3*mu2s);
if X2 <= 0
    ss.washout = true;
    ss = fill_washout(ss, S1in, S2in);
    return
end

ss.S1  = S1;
ss.X1  = X1;
ss.S2  = S2;
ss.X2  = X2;
ss.mu1 = mu1s;
ss.mu2 = mu2s;
ss.fT  = fT;
ss.qM  = a.k6 * mu2s * X2;        % mmol CH4/(L d)
end

% =====================================================================
function ss = fill_washout(ss, S1in, S2in)
ss.S1 = S1in;  ss.X1 = 0;
ss.S2 = S2in;  ss.X2 = 0;
ss.mu1 = 0;    ss.mu2 = 0;  ss.fT = 0;  ss.qM = 0;
end
