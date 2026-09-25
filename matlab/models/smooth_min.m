function m = smooth_min(a, b, eps_s)
%SMOOTH_MIN Differentiable minimum of two quantities.
%
%   M = SMOOTH_MIN(A, B, EPS_S) returns
%       0.5*(A + B - sqrt((A-B)^2 + EPS_S^2))
%   which approaches min(A,B) as EPS_S -> 0 but, unlike MIN, has a
%   continuous derivative where the two arguments cross.
%
%   Used wherever a supply saturates a demand -- available heat against
%   heat demand, stored feedstock against feed demand. The crossing point
%   is exactly where such a system spends most of its time, so a kink there
%   is a kink the solver will keep hitting.
%
%   See also SMOOTH_STEP, DIGESTER_THERMOSTAT.

if nargin < 3 || isempty(eps_s)
    eps_s = 1e-3;
end
m = 0.5*(a + b - sqrt((a - b).^2 + eps_s^2));
end
