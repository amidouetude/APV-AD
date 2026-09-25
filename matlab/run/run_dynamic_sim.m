function R = run_dynamic_sim(clim, site, design, P, D, opts)
%RUN_DYNAMIC_SIM Integrate the coupled APV-AD model over one or more years.
%
%   R = RUN_DYNAMIC_SIM(CLIM, SITE, DESIGN, P, D) integrates the model with
%   a one-year spin-up and returns the final year.
%
%   R = RUN_DYNAMIC_SIM(..., OPTS) accepts a struct with fields
%     mode        'apv' (default) or 'reference' (open-field, unshaded)
%     n_years     total years integrated, including spin-up (default 2)
%     out_dt_d    output sampling interval in days (default 1/24)
%     solver      'ode15s' (default), 'ode23t', 'ode23tb'
%     verbose     print progress (default true)
%     x0          override the initial state
%
%   Output R contains
%     R.t        n x 1   time in days within the reported year
%     R.time     n x 1   datetime
%     R.X        n x m   states on the output grid
%     R.out      struct of algebraic outputs, one field per quantity
%     R.annual   struct of annual totals and averages
%     R.F        the forcing struct used
%     R.sol      the raw solver solution structure
%     R.stats    solver statistics and wall-clock time
%
%   Why a stiff solver
%   ------------------
%   The state vector spans time constants from about 12 minutes (PV module
%   thermal mass) through a few days (digester temperature, methanogen
%   growth) to a full season (thermal time, feedstock store). That is four
%   orders of magnitude in one system, which is the definition of a stiff
%   problem: an explicit method would be forced to the smallest time
%   constant for the whole year purely for stability, not accuracy. ODE15S
%   is variable-order NDF and steps on accuracy alone.
%
%   MaxStep is capped at the forcing interval so the solver cannot step
%   over a day of weather -- the interpolants are smooth enough that it
%   would otherwise happily take multi-day steps through the night.
%
%   See also APVAD_RHS, POSTPROCESS_SIM, BUILD_FORCING.

if nargin < 6, opts = struct(); end
opts = set_default(opts, 'mode',     'apv');
opts = set_default(opts, 'n_years',  2);
opts = set_default(opts, 'out_dt_d', 1/24);
opts = set_default(opts, 'solver',   'ode15s');
opts = set_default(opts, 'verbose',  true);

% ---------------------------------------------------------------------
% Forcing and initial condition
% ---------------------------------------------------------------------
F = build_forcing(clim, site, design, P, D, opts.mode);

if isfield(opts, 'x0') && ~isempty(opts.x0)
    x0 = opts.x0(:);
else
    x0 = apvad_initial_state(clim, F, D, P);
end

tEnd  = opts.n_years * F.t_period;
tspan = [0, tEnd];

ixs = apvad_states();

odeOpts = odeset( ...
    'RelTol',  D.num.rel_tol, ...
    'AbsTol',  D.num.abs_tol, ...
    'MaxStep', D.num.max_step_d, ...
    'JPattern', ixs.JPattern, ...
    'Stats',   'off');

if isfield(D.num,'nonnegative') && D.num.nonnegative
    odeOpts = odeset(odeOpts, 'NonNegative', nonneg_indices());
end

if opts.verbose
    fprintf(['[run_dynamic_sim] %s | %s | %d year(s), spin-up %d | ' ...
             'solver %s\n'], site.name, opts.mode, opts.n_years, ...
             opts.n_years-1, opts.solver);
end

rhs = @(t, x) apvad_rhs(t, x, F, D, P);

tic;
switch lower(opts.solver)
    case 'ode15s',  sol = ode15s(rhs, tspan, x0, odeOpts);
    case 'ode23t',  sol = ode23t(rhs, tspan, x0, odeOpts);
    case 'ode23tb', sol = ode23tb(rhs, tspan, x0, odeOpts);
    otherwise
        error('run_dynamic_sim:unknownSolver', ...
              'Unsupported solver "%s".', opts.solver);
end
wall = toc;

% ---------------------------------------------------------------------
% Sample the reported (final) year
% ---------------------------------------------------------------------
t0_report = (opts.n_years - 1) * F.t_period;
t_grid    = (t0_report : opts.out_dt_d : tEnd)';
if t_grid(end) < tEnd, t_grid(end+1) = tEnd; end

X = deval(sol, t_grid)';

R.t     = t_grid - t0_report;
R.t_abs = t_grid;
R.X     = X;
R.F     = F;
R.sol   = sol;
R.site  = site;
R.design = design;
R.mode  = opts.mode;

yr = year(clim.time(1));
R.time = datetime(yr,1,1) + days(R.t);

R.stats.wall_s     = wall;
R.stats.n_steps    = numel(sol.x);
R.stats.n_years    = opts.n_years;
R.stats.solver     = opts.solver;
R.stats.mean_step_d = mean(diff(sol.x));

if opts.verbose
    fprintf(['[run_dynamic_sim] done in %.1f s, %d solver steps ' ...
             '(mean step %.4f d)\n'], wall, R.stats.n_steps, ...
             R.stats.mean_step_d);
end

% ---------------------------------------------------------------------
% Algebraic outputs and annual aggregates
% ---------------------------------------------------------------------
R = postprocess_sim(R, D, P);
end

% =====================================================================
function s = set_default(s, f, v)
if ~isfield(s, f) || isempty(s.(f))
    s.(f) = v;
end
end

% =====================================================================
function idx = nonneg_indices()
%NONNEG_INDICES States that must never go negative.
%   Concentrations, stores and accumulators are constrained; temperatures
%   and the depletion Dr are not, because Dr < 0 legitimately represents
%   water above field capacity awaiting percolation and temperatures can be
%   below zero at Freiburg and Konya in winter.
ix = apvad_states();
idx = [ix.S1, ix.X1, ix.S2, ix.X2, ix.Z, ix.C, ...
       ix.GDD, ix.B_dm, ix.M_res, ix.accumulators];
end
