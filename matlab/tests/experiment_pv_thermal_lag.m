function R = experiment_pv_thermal_lag(siteName, nDays)
%EXPERIMENT_PV_THERMAL_LAG What the PV module's thermal memory is worth.
%
%   R = EXPERIMENT_PV_THERMAL_LAG() runs 30 days at Konya twice -- once
%   with the algebraic Faiman module temperature and once with the
%   first-order thermal lag -- and reports the difference in PV yield,
%   module temperature and runtime.
%
%   R = EXPERIMENT_PV_THERMAL_LAG(SITE, NDAYS) for another site or horizon.
%
%   Why this experiment exists
%   --------------------------
%   D.pv.dynamic defaults to false, and a default that switches off a
%   feature deserves a measurement rather than an assertion. The module
%   time constant is about twelve minutes, two orders of magnitude faster
%   than the digester's, so resolving it drags the stiff solver down to
%   roughly one-minute steps for the entire simulated year.
%
%   The question is what that buys. Against hourly forcing, very little for
%   an annual energy total: the module spends the year alternately lagging
%   a rising sun and a falling one, and the two errors very nearly cancel
%   in the integral. It buys a great deal for anything that looks at the
%   module's response in time -- sub-hourly yield, ramp rates, thermal
%   cycling -- which is why the capability stays in the model rather than
%   being removed.
%
%   Report R.dPV_pct and R.T_rms_K in the manuscript alongside the claim
%   that the algebraic form is used for annual runs.
%
%   See also PV_MODULE, RUN_DYNAMIC_SIM, APVAD_DYN_PARAMS.

if nargin < 1 || isempty(siteName), siteName = 'Konya'; end
if nargin < 2 || isempty(nDays),    nDays = 30;         end

thisDir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(thisDir));
projRoot = fileparts(thisDir);

P    = apvad_params();
D    = apvad_dyn_params();
site = apvad_sites(siteName);
ix   = apvad_states();

w = warning('off','load_climate:noPrecipitation');
clim = load_climate(siteName, fullfile(projRoot,'outputs','csv'));
warning(w);

design = struct('beta_deg',25,'d_row_m',6,'H_m_m',2.5, ...
                'V_dig_m3',10,'HRT_days',30,'f_PV_heat',0.15);

% Start in high summer, where irradiance is strongest and the lag largest.
t0 = 170;
tspan = [t0, t0 + nDays];
grid_t = (tspan(1) : 1/60 : tspan(2))';       % one-minute output grid

out = struct();
for mode = ["static", "dynamic"]
    Dm = D;
    Dm.pv.dynamic = strcmp(mode, "dynamic");

    F  = build_forcing(clim, site, design, P, Dm, 'apv');
    x0 = apvad_initial_state(clim, F, Dm, P);
    x0(ix.T_cell) = F.T_amb(t0);

    opt = odeset('RelTol',Dm.num.rel_tol, 'AbsTol',Dm.num.abs_tol, ...
                 'MaxStep',Dm.num.max_step_d, 'JPattern',ix.JPattern);

    tic;
    sol = ode15s(@(t,x) apvad_rhs(t,x,F,Dm,P), tspan, x0, opt);
    wall = toc;

    X = deval(sol, grid_t)';
    n = numel(grid_t);
    Tc = zeros(n,1); Pw = zeros(n,1);
    for i = 1:n
        [~, o] = apvad_rhs(grid_t(i), X(i,:)', F, Dm, P);
        Tc(i) = o.T_cell;
        Pw(i) = o.P_array_W;
    end

    out.(mode).E_kWh   = X(end, ix.E_pv) - X(1, ix.E_pv);
    out.(mode).T_cell  = Tc;
    out.(mode).P_W     = Pw;
    out.(mode).wall_s  = wall;
    out.(mode).n_steps = numel(sol.x);
end

R.t          = grid_t;
R.static     = out.static;
R.dynamic    = out.dynamic;
R.dPV_pct    = 100*(out.dynamic.E_kWh - out.static.E_kWh) ...
               / max(out.static.E_kWh, eps);
R.T_rms_K    = sqrt(mean((out.dynamic.T_cell - out.static.T_cell).^2));
R.T_max_dev_K= max(abs(out.dynamic.T_cell - out.static.T_cell));
R.speedup    = out.dynamic.wall_s / max(out.static.wall_s, eps);
R.site       = siteName;
R.n_days     = nDays;

fprintf('\n=== PV thermal lag, %s, %d days from day %d ===\n', ...
        siteName, nDays, t0);
fprintf('  PV energy, algebraic Faiman : %10.1f kWh/ha\n', out.static.E_kWh);
fprintf('  PV energy, thermal lag      : %10.1f kWh/ha\n', out.dynamic.E_kWh);
fprintf('  difference                  : %10.3f %%\n', R.dPV_pct);
fprintf('  module temperature RMS diff : %10.2f K\n', R.T_rms_K);
fprintf('  largest instantaneous diff  : %10.2f K\n', R.T_max_dev_K);
fprintf('  solver steps  %d -> %d\n', out.static.n_steps, out.dynamic.n_steps);
fprintf('  runtime       %.1f s -> %.1f s  (x%.1f)\n\n', ...
        out.static.wall_s, out.dynamic.wall_s, R.speedup);
end
