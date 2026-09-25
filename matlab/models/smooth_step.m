function s = smooth_step(z, w)
%SMOOTH_STEP Differentiable approximation of the Heaviside step.
%
%   S = SMOOTH_STEP(Z, W) returns 0.5*(1 + tanh(Z/W)), which goes from 0 to
%   1 over a transition of width roughly W around Z = 0.
%
%   Every switching decision in this model -- start irrigating, stop
%   drawing feed from an empty store, turn the digester heater on -- is
%   expressed with this function rather than an IF statement. A variable-
%   step stiff solver estimates the local error from a smooth Taylor
%   expansion of the right-hand side; a genuine discontinuity makes that
%   estimate meaningless, so the step size collapses and the integration
%   either crawls or silently steps straight over the event. The
%   alternatives are event functions with restarts (correct, but heavy for
%   switches that fire hundreds of times a year) or smoothing. W is chosen
%   narrow enough to be physically negligible and wide enough for the
%   solver to resolve; see D.num.switch_width.
%
%   See also SOIL_WATER, FEED_SCHEDULER, DIGESTER_THERMOSTAT.

if nargin < 2 || isempty(w) || w <= 0
    w = 1e-3;
end
s = 0.5*(1 + tanh(z./w));
end
