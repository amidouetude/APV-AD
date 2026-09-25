function figs = plot_apvad_dynamics(OUT, figDir, zoomStartDay)
%PLOT_APVAD_DYNAMICS Publication figures for the dynamic APV-AD model.
%
%   FIGS = PLOT_APVAD_DYNAMICS(OUT) draws two figures from the struct
%   returned by MAIN_APVAD_DYNAMIC and returns their handles.
%
%   FIGS = PLOT_APVAD_DYNAMICS(OUT, FIGDIR) also writes 300 dpi PNGs.
%   FIGS = PLOT_APVAD_DYNAMICS(OUT, FIGDIR, ZOOMSTARTDAY) sets the start
%   day of the ten-day detail window (default: the coldest ten days of the
%   year, where the heating system is worked hardest).
%
%   Figure 1 -- annual dynamics, eight panels of daily means.
%   Figure 2 -- a ten-day window at hourly resolution, showing a heat
%   supply that only exists in daylight against a demand that does not stop
%   at sunset.
%
%   Chart conventions
%   -----------------
%   ONE MEASUREMENT SCALE PER PANEL. No panel carries a secondary y-axis,
%   and no panel plots series whose magnitudes differ by more than about an
%   order of magnitude: PV generation runs to 200 kW while the digester
%   asks for well under one, so putting them together would draw the heat
%   demand as a flat line on the axis and invite the reader to conclude it
%   is zero. They get their own panels instead.
%
%   Y-AXES INCLUDE THEIR MEANINGFUL ZERO OR REFERENCE BAND. Several states
%   here are, correctly, almost constant -- a well-fed digester held at
%   setpoint produces very nearly constant methane. On an auto-scaled axis
%   a one-percent ripple fills the panel and reads as instability. Each
%   panel is therefore scaled against something physical: zero for rates
%   and concentrations, the inhibition threshold for VFA, the operational
%   band for pH, total available water for depletion.
%
%   Series colours are a fixed three-hue categorical set (blue, orange,
%   aqua) validated for colour-vision deficiency across all pairs, assigned
%   in fixed order and never cycled. Reference levels are recessive grey
%   dashed rules, since they are annotations rather than measurements.
%   Every panel with two or more series carries a legend.
%
%   Annual panels show DAILY MEANS: 8760 hourly points per series on a
%   150 mm page is ink, not information. Figure 2 supplies hourly detail
%   where it carries the argument.
%
%   See also MAIN_APVAD_DYNAMIC, PLOT_FROM_SAVED, POSTPROCESS_SIM.

if nargin < 2, figDir = ''; end

A   = OUT.apv;
Rf  = OUT.ref;
o   = A.out;
D   = OUT.D;
P   = OUT.P;
ann = A.annual;

% ---- Validated categorical palette (fixed order, never cycled) -------
c1 = [0.165 0.471 0.839];    % #2a78d6 blue
c2 = [0.922 0.408 0.204];    % #eb6834 orange
c3 = [0.106 0.686 0.478];    % #1baf7a aqua
cg = [0.55  0.55  0.55];     % recessive grey for reference rules
LW = 1.6;
YEAR = [0 365];

t = A.t;                     % days within the reported year

% =====================================================================
% Figure 1 -- annual dynamics (daily means)
% =====================================================================
figs(1) = figure('Name', sprintf('APV-AD annual dynamics -- %s', ann.site), ...
                 'Color','w', 'Position',[80 40 1180 940]);
T = tiledlayout(figs(1), 4, 2, 'TileSpacing','compact', 'Padding','compact');
title(T, sprintf(['Dynamic APV-AD model -- %s   (daily means, ' ...
                  'periodic steady state)'], ann.site), ...
      'FontWeight','bold','FontSize',13);

[td, Tdig] = dayavg(t, o.T_dig);
[~,  Tamb] = dayavg(t, o.T_amb);

% --- A. Digester temperature ---
ax = nexttile; hold(ax,'on'); grid(ax,'on');
plot(ax, td, Tdig, 'Color',c1, 'LineWidth',LW, 'DisplayName','Digester');
plot(ax, td, Tamb, 'Color',c2, 'LineWidth',LW, 'DisplayName','Ambient air');
yline(ax, D.thermal.T_set_degC, '--', 'Color',cg, 'LineWidth',1, ...
      'DisplayName','Setpoint 37 \circC');
yline(ax, P.T_dig_min_degC, ':', 'Color',cg, 'LineWidth',1, ...
      'DisplayName',sprintf('Minimum %.0f \\circC', P.T_dig_min_degC));
style(ax, 'A  Digester temperature', 'Temperature (\circC)', YEAR);
legend(ax,'Location','southeast','Box','off','FontSize',8);
note(ax, sprintf('%.0f %% of the year above 35 \\circC', ...
                 100*ann.frac_time_mesophilic));

% --- B. Methane production ---
ax = nexttile; hold(ax,'on'); grid(ax,'on');
[~, qch4] = dayavg(t, o.Q_CH4_Nm3_d);
plot(ax, td, qch4, 'Color',c1, 'LineWidth',LW);
style(ax, 'B  Methane production', 'Q_{CH_4} (Nm^3 d^{-1})', YEAR);
ylim(ax, [0, max(qch4)*1.3]);
note(ax, sprintf('annual %.0f Nm^3 ha^{-1}', ann.CH4_Nm3_ha));

% --- C. Volatile fatty acids ---
ax = nexttile; hold(ax,'on'); grid(ax,'on');
[~, s2] = dayavg(t, o.S2);
plot(ax, td, s2, 'Color',c1, 'LineWidth',LW, 'DisplayName','VFA S_2');
yline(ax, ann.S2_crit_mmolL, '--', 'Color',cg, 'LineWidth',1, ...
      'DisplayName',sprintf('Haldane peak %.0f mmol L^{-1}', ann.S2_crit_mmolL));
style(ax, 'C  Volatile fatty acids', 'S_2 (mmol L^{-1})', YEAR);
ylim(ax, [0, ann.S2_crit_mmolL*1.25]);
legend(ax,'Location','northeast','Box','off','FontSize',8);

% --- D. Digester pH ---
ax = nexttile; hold(ax,'on'); grid(ax,'on');
[~, ph] = dayavg(t, o.pH);
plot(ax, td, ph, 'Color',c1, 'LineWidth',LW, 'DisplayName','pH');
yline(ax, 6.8, '--', 'Color',cg, 'LineWidth',1, ...
      'DisplayName','Souring threshold');
style(ax, 'D  Digester pH', 'pH (-)', YEAR);
ylim(ax, [6.4, 8.0]);
legend(ax,'Location','northeast','Box','off','FontSize',8);

% --- E. PV generation ---
% Exported power is within a fraction of a percent of generation at this
% design point, and the slice diverted to the digester is three orders of
% magnitude smaller; both would be invisible here. Panel F carries the
% heat side on its own scale.
ax = nexttile; hold(ax,'on'); grid(ax,'on');
[~, ppv] = dayavg(t, o.P_array_W/1000);
plot(ax, td, ppv, 'Color',c1, 'LineWidth',LW);
style(ax, 'E  PV array generation', 'Power (kW ha^{-1})', YEAR);
ylim(ax, [0, max(ppv)*1.15]);
note(ax, sprintf('annual %.0f MWh ha^{-1}; %.2f %% diverted to heat', ...
                 ann.PV_total_MWh_ha, ...
                 100*ann.heat_from_pv_kWh_ha/max(ann.PV_total_kWh_ha,eps)));

% --- F. Digester heat balance ---
ax = nexttile; hold(ax,'on'); grid(ax,'on');
[~, qd] = dayavg(t, o.Q_demand_W/1000);
[~, qp] = dayavg(t, o.Q_pv_W/1000);
[~, qb] = dayavg(t, o.Q_boiler_W/1000);
plot(ax, td, qd, 'Color',c1, 'LineWidth',LW, 'DisplayName','Demand');
plot(ax, td, qp, 'Color',c2, 'LineWidth',LW, 'DisplayName','From PV');
plot(ax, td, qb, 'Color',c3, 'LineWidth',LW, 'DisplayName','From biogas boiler');
style(ax, 'F  Digester heating', 'Heat rate (kW)', YEAR);
ylim(ax, [0, max(qd)*1.35]);
legend(ax,'Location','northeast','Box','off','FontSize',8);
note(ax, sprintf('demand covered %.0f %%', 100*ann.heat_coverage));

% --- G. Root-zone water depletion ---
ax = nexttile; hold(ax,'on'); grid(ax,'on');
[~, dr] = dayavg(t, o.Dr);
plot(ax, td, dr, 'Color',c1, 'LineWidth',LW, 'DisplayName','Depletion D_r');
yline(ax, ann.RAW_mm, '--', 'Color',cg, 'LineWidth',1, ...
      'DisplayName','RAW (stress onset)');
yline(ax, ann.TAW_mm, ':', 'Color',cg, 'LineWidth',1, 'DisplayName','TAW');
style(ax, 'G  Root-zone water depletion', 'D_r (mm)', YEAR);
ylim(ax, [0, ann.TAW_mm*1.08]);
legend(ax,'Location','northeast','Box','off','FontSize',8);
xlabel(ax, 'Day of year');

% --- H. Crop biomass, agrivoltaic vs open field ---
ax = nexttile; hold(ax,'on'); grid(ax,'on');
[~, bapv] = dayavg(t, o.B_dm/1000);
[~, bref] = dayavg(Rf.t, Rf.out.B_dm/1000);
plot(ax, td, bapv, 'Color',c1, 'LineWidth',LW, 'DisplayName','Under array');
plot(ax, td, bref, 'Color',c2, 'LineWidth',LW, 'DisplayName','Open field');
style(ax, 'H  Standing crop biomass', 'Dry matter (t ha^{-1})', YEAR);
ylim(ax, [0, max(bref)*1.25]);
legend(ax,'Location','northwest','Box','off','FontSize',8);
xlabel(ax, 'Day of year');

% =====================================================================
% Figure 2 -- ten-day detail window
% =====================================================================
if nargin < 3 || isempty(zoomStartDay)
    zoomStartDay = coldest_window(t, o);
end
w  = (t >= zoomStartDay) & (t <= zoomStartDay + 10);
XW = [zoomStartDay, zoomStartDay + 10];

figs(2) = figure('Name', sprintf('APV-AD detail window -- %s', ann.site), ...
                 'Color','w', 'Position',[120 60 1060 900]);
T2 = tiledlayout(figs(2), 4, 1, 'TileSpacing','compact','Padding','compact');
title(T2, sprintf(['%s -- ten days from day %d, hourly: a heat supply ' ...
                   'that only exists in daylight'], ...
                   ann.site, round(zoomStartDay)), ...
      'FontWeight','bold','FontSize',12);

ax = nexttile; hold(ax,'on'); grid(ax,'on');
plot(ax, t(w), o.T_dig(w), 'Color',c1, 'LineWidth',LW, 'DisplayName','Digester');
plot(ax, t(w), o.T_amb(w), 'Color',c2, 'LineWidth',LW, 'DisplayName','Ambient air');
yline(ax, D.thermal.T_set_degC, '--','Color',cg,'LineWidth',1, ...
      'DisplayName','Setpoint');
style(ax, 'Digester and ambient temperature', 'Temperature (\circC)', XW);
legend(ax,'Location','east','Box','off','FontSize',8);

% Heat demand and PV generation differ by more than two orders of
% magnitude, so they occupy separate panels rather than one axis on which
% the digester's demand would be indistinguishable from zero.
ax = nexttile; hold(ax,'on'); grid(ax,'on');
% Demand is drawn as a wide, pale envelope beneath the two supply traces.
% When coverage is complete the supplies sum to it exactly and would
% otherwise hide it entirely, leaving the reader unable to see that the
% handover between PV and boiler is a handover rather than a shortfall.
plot(ax, t(w), o.Q_demand_W(w), 'Color',[c1 0.35], 'LineWidth',4*LW, ...
     'DisplayName','Total demand');
plot(ax, t(w), o.Q_pv_W(w),     'Color',c2,'LineWidth',LW,'DisplayName','Supplied by PV');
plot(ax, t(w), o.Q_boiler_W(w), 'Color',c3,'LineWidth',LW,'DisplayName','Supplied by boiler');
style(ax, 'Digester heat: demand and where it comes from', 'Heat rate (W)', XW);
ylim(ax, [0, max(o.Q_demand_W(w))*1.35]);
legend(ax,'Location','east','Box','off','FontSize',8);

ax = nexttile; hold(ax,'on'); grid(ax,'on');
plot(ax, t(w), o.P_array_W(w)/1000, 'Color',c1,'LineWidth',LW);
style(ax, 'PV array generation', 'Power (kW ha^{-1})', XW);
ylim(ax, [0, max(o.P_array_W(w)/1000)*1.15]);
note(ax, 'zero every night; the digester''s demand is not');

ax = nexttile; hold(ax,'on'); grid(ax,'on');
plot(ax, t(w), o.Q_CH4_Nm3_d(w), 'Color',c1,'LineWidth',LW);
style(ax, 'Methane production rate', 'Q_{CH_4} (Nm^3 d^{-1})', XW);
ylim(ax, [0, max(o.Q_CH4_Nm3_d(w))*1.3]);
xlabel(ax, 'Day of year');

% =====================================================================
% Save
% =====================================================================
if ~isempty(figDir)
    if ~isfolder(figDir), mkdir(figDir); end
    f1 = fullfile(figDir, sprintf('dyn_annual_%s.png', ann.site));
    f2 = fullfile(figDir, sprintf('dyn_detail_%s.png', ann.site));
    exportgraphics(figs(1), f1, 'Resolution', 300);
    exportgraphics(figs(2), f2, 'Resolution', 300);
    fprintf('saved: %s\nsaved: %s\n', f1, f2);
end
end

% =====================================================================
function style(ax, ttl, ylab, xl)
title(ax, ttl, 'FontWeight','normal', 'FontSize',10);
ylabel(ax, ylab, 'FontSize',9);
ax.FontSize   = 8;
ax.Box        = 'off';
ax.GridColor  = [0.8 0.8 0.8];
ax.GridAlpha  = 0.6;
ax.XLim       = xl;
end

% =====================================================================
function note(ax, txt)
%NOTE Small recessive annotation in the upper-left of a panel.
text(ax, 0.015, 0.93, txt, 'Units','normalized', ...
     'FontSize',8, 'Color',[0.35 0.35 0.35], ...
     'VerticalAlignment','top', 'Interpreter','tex');
end

% =====================================================================
function [td, vd] = dayavg(t, v)
%DAYAVG Daily means of an hourly series.
d = floor(t);
[u, ~, g] = unique(d);
vd = accumarray(g, v, [], @mean);
td = u;
end

% =====================================================================
function d0 = coldest_window(t, o)
%COLDEST_WINDOW Start day of the ten-day stretch with the lowest mean
%   ambient temperature -- where the heating system is worked hardest and
%   the coupling is most visible.
d = floor(t);
[u, ~, g] = unique(d);
Td = accumarray(g, o.T_amb, [], @mean);
m  = movmean(Td, 10);
m(u > max(u) - 11) = inf;
[~, i] = min(m);
d0 = u(i);
end
