% smoke_test.m -- fast structural check of the dynamic APV-AD model.
% Run this before the full simulation: it exercises every sub-model on a
% short horizon so that syntax errors, unit slips and shape mismatches
% surface in seconds rather than after a full-year integration.

clear; clc;
thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(thisDir)));
projRoot = fileparts(fileparts(thisDir));

fprintf('--- 1. configuration ---\n');
P = apvad_params();
D = apvad_dyn_params();
ix = apvad_states();
site = apvad_sites('Konya');
fprintf('states: %d | site: %s (lat %.2f)\n', ix.n, site.name, site.lat);

fprintf('--- 2. climate ---\n');
clim = load_climate('Konya', fullfile(projRoot,'outputs','csv'));
fprintf('rows: %d | dt: %.2f h | t: %.3f .. %.3f d | season hours: %d\n', ...
    numel(clim.t_day), clim.dt_h, clim.t_day(1), clim.t_day(end), ...
    sum(clim.in_season));
fprintf('GHI annual: %.0f kWh/m2 | Tmean: %.2f degC\n', ...
    sum(clim.GHI)/1000, mean(clim.T2m));

fprintf('--- 3. forcing ---\n');
design = struct('beta_deg',20.71,'d_row_m',3.030,'H_m_m',1.501, ...
                'V_dig_m3',4.317,'HRT_days',22.51,'f_PV_heat',0.0513);
F = build_forcing(clim, site, design, P, D, 'apv');
fprintf('n_modules: %d | GCR: %.3f | Q_in: %.4f m3/d | A_wall: %.2f m2\n', ...
    F.n_modules, F.GCR, F.Q_in_m3_d, F.A_wall_m2);
fprintf('c_VS_in: %.2f g/L | nominal OLR: %.3f kgVS/m3/d | VS_ann: %.0f kg\n', ...
    F.c_VS_in_gL, F.c_VS_in_gL/design.HRT_days, F.VS_ann_nom_kg);
fprintf('G_eff annual: %.0f kWh/m2 | F_shad mean: %.3f\n', ...
    trapz(F.t_day, F.G_eff(F.t_day))*24/1000, mean(F.F_shad(F.t_day)));
fprintf('heater: design %.0f W | Q_max %.0f W | Kp %.0f W/K | T_amb_design %.1f degC\n', ...
    F.Q_design_W, F.Q_max_W, F.Kp_W_per_K, F.T_amb_design_degC);

fprintf('--- 4. AM2 steady state at 37 degC ---\n');
S1in = D.feed.f_S1*D.feed.COD_per_VS*F.c_VS_in_gL;
ss = am2_steady_state(1/design.HRT_days, S1in, D.feed.S2in_mmolL, 37, D, P);
fprintf('washout: %d | S1*: %.3f | X1*: %.3f | S2*: %.2f | X2*: %.3f\n', ...
    ss.washout, ss.S1, ss.X1, ss.S2, ss.X2);
fprintf('qM: %.2f mmol/(L d)\n', ss.qM);

fprintf('--- 5. BMP calibration ---\n');
[Dc, ci] = calibrate_am2_to_bmp(D, P, design.HRT_days, F.c_VS_in_gL);
fprintf('BMP model %.1f -> target %.1f (k6 %.1f -> %.1f, scale %.4f)\n', ...
    ci.BMP_model, ci.BMP_target, ci.k6_before, ci.k6_after, ci.scale);

fprintf('--- 6. single RHS evaluation ---\n');
x0 = apvad_initial_state(clim, F, Dc, P);
[dx, o] = apvad_rhs(100, x0, F, Dc, P);
fprintf('dx finite: %d | numel: %d\n', all(isfinite(dx)), numel(dx));
fprintf('T_cell_ss %.2f | P_array %.1f W | Q_CH4 %.4f Nm3/d | pH %.2f\n', ...
    o.T_cell_ss, o.P_array_W, o.Q_CH4_Nm3_d, o.pH);
fprintf('tau_dig %.2f d | Q_loss %.1f W | Q_demand %.1f W\n', ...
    o.tau_dig_d, o.Q_loss_W, o.Q_demand_W);

fprintf('--- 7. 20-day integration and timing probe ---\n');
opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1/24, ...
              'JPattern',ix.JPattern);
tic;
sol = ode15s(@(t,x) apvad_rhs(t,x,F,Dc,P), [0 20], x0, opts);
el = toc;
fprintf('ok: %d steps in %.2f s (mean step %.4f d)\n', ...
    numel(sol.x), el, mean(diff(sol.x)));
fprintf('extrapolated one-year runtime: %.0f s\n', el*365/20);
xe = sol.y(:,end);
fprintf('T_dig %.2f | S1 %.3f | S2 %.2f | X2 %.3f | E_pv %.1f kWh\n', ...
    xe(ix.T_dig), xe(ix.S1), xe(ix.S2), xe(ix.X2), xe(ix.E_pv));

fprintf('\nSMOKE TEST COMPLETE\n');
