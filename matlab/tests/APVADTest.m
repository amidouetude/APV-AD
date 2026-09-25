classdef APVADTest < matlab.unittest.TestCase
%APVADTEST Unit tests for the dynamic APV-AD model.
%
%   Run with
%       results = runtests('APVADTest')
%   or, from the matlab/ folder,
%       run_all_tests
%
%   The tests are chosen to catch the failure modes that a plausible-looking
%   annual number would hide. Several check the model against a CLOSED FORM
%   rather than against a stored expected value, so they stay meaningful if
%   parameters change: the angle of incidence against a trigonometric
%   identity at solar noon, the PV thermal ODE against the Faiman
%   expression it must relax to, and the integrated AM2 digester against
%   its own analytical equilibrium.
%
%   See also RUN_ALL_TESTS, SMOKE_TEST.

    properties
        P
        D
        site
        clim
        design
    end

    methods (TestClassSetup)
        function setupPaths(tc)
            here = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(here));
            tc.P = apvad_params();
            tc.D = apvad_dyn_params();
            tc.site = apvad_sites('Konya');
            tc.design = struct('beta_deg',20.71,'d_row_m',3.030, ...
                               'H_m_m',1.501,'V_dig_m3',4.317, ...
                               'HRT_days',22.51,'f_PV_heat',0.0513);
        end
    end

    methods (Access = private)
        function c = getClimate(tc)
            if isempty(tc.clim)
                projRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
                w = warning('off','load_climate:noPrecipitation');
                tc.clim = load_climate('Konya', ...
                            fullfile(projRoot,'outputs','csv'));
                warning(w);
            end
            c = tc.clim;
        end
    end

% =====================================================================
    methods (Test)  % --- forcing layer -------------------------------

    function testClimateAxisIsStrictlyIncreasing(tc)
    %   The TMY file splices months from different calendar years, so its
    %   raw timestamps run backwards at month boundaries. The loader must
    %   rebuild a single monotone year; if it does not, griddedInterpolant
    %   refuses the data and nothing downstream runs.
        c = tc.getClimate();
        tc.verifyTrue(all(diff(c.t_day) > 0), ...
            'Climate time axis must be strictly increasing.');
        tc.verifyEqual(numel(c.t_day), 8760, ...
            'A typical meteorological year must hold 8760 hourly rows.');
        tc.verifyLessThan(max(c.t_day), 365);
    end

    function testSolarNoonElevationMatchesIdentity(tc)
    %   At true solar noon the elevation is exactly 90 - |lat - dec|.
        doy = (1:20:361)';
        sun = solar_angles(doy, zeros(size(doy)), tc.site.lat, tc.site.lon);
        % Re-evaluate at the clock hour that puts the hour angle at zero.
        hNoon = 12 - tc.site.lon/15 - sun.eot_min/60;
        sunN  = solar_angles(doy, hNoon, tc.site.lat, tc.site.lon);
        tc.verifyEqual(sunN.hour_angle_deg, zeros(size(doy)), 'AbsTol', 1e-9);
        expected = 90 - abs(tc.site.lat - sunN.declination_deg);
        tc.verifyEqual(sunN.elevation_deg, expected, 'AbsTol', 1e-8);
    end

    function testAngleOfIncidenceIdentityAtNoon(tc)
    %   With the hour angle at zero the Duffie-Beckman expression collapses
    %   to cos(theta_i) = cos(lat - beta - dec). This is the identity the
    %   original Python code violated by assuming the sun stays due south
    %   all day, so it is worth pinning down.
        dec  = (-23:5:23)';
        beta = 25;
        lat  = tc.site.lat;
        poa = poa_irradiance(zeros(size(dec)), ones(size(dec)), ...
                             zeros(size(dec)), dec, zeros(size(dec)), ...
                             zeros(size(dec)), beta, lat, tc.P);
        expected = max(cosd(lat - beta - dec), 0);
        tc.verifyEqual(poa.cos_theta_i, expected, 'AbsTol', 1e-10);
    end

    function testShadingFractionBounds(tc)
        elev = (-10:1:90)';
        [Fs, GCR] = shading_fraction(elev, 25, 4, 2, tc.P);
        tc.verifyTrue(all(Fs >= 0 & Fs <= 1));
        tc.verifyTrue(all(Fs(elev <= 0) == 0), ...
            'There is no ground shading when the sun is below the horizon.');
        tc.verifyGreaterThan(GCR, 0);
        % Shading falls as the sun rises.
        hi = Fs(elev > 5);
        tc.verifyTrue(all(diff(hi) <= 1e-12));
    end

    end

% =====================================================================
    methods (Test)  % --- photovoltaics -------------------------------

    function testPVThermalReducesToFaiman(tc)
    %   The dynamic module temperature must have the Faiman expression as
    %   its exact steady state. If this drifts, the dynamic and static PV
    %   yields are no longer comparable and the whole benchmark is
    %   confounded.
        Dd = tc.D;  Dd.pv.dynamic = true;  Dd.pv.include_Pelec = false;
        G = 800; Ta = 20; WS = 2;
        Tss = Ta + G/(tc.P.U0 + tc.P.U1*WS);

        pv = pv_module(G, Ta, WS, Tss, 1, Dd, tc.P);
        tc.verifyEqual(pv.dTcell_dt, 0, 'AbsTol', 1e-10, ...
            'At the Faiman temperature the thermal derivative must vanish.');
        tc.verifyEqual(pv.T_cell_ss, Tss, 'AbsTol', 1e-12);

        % Relaxing from a cold start must converge to the same value.
        % The integration tolerance has to be tighter than the agreement
        % being asserted, otherwise the test measures ODE45's default
        % RelTol of 1e-3 rather than the model.
        f = @(t,T) pv_thermal_rhs(T, G, Ta, WS, Dd, tc.P);
        s = ode45(f, [0 0.5], Ta, odeset('RelTol',1e-10,'AbsTol',1e-12));
        tc.verifyEqual(s.y(end), Tss, 'RelTol', 1e-6);
    end

    function testPVPowerTemperatureDerating(tc)
    %   Power must fall with module temperature at the stated coefficient.
        G = 1000; Ta = 25; WS = 1;
        p25 = pv_module(G, Ta, WS, 25, 1, tc.D, tc.P);
        p45 = pv_module(G, Ta, WS, 45, 1, tc.D, tc.P);
        ratio = p45.P_module_W / p25.P_module_W;
        tc.verifyEqual(ratio, 1 + tc.P.beta_p*20, 'RelTol', 1e-12);
    end

    end

% =====================================================================
    methods (Test)  % --- anaerobic digestion -------------------------

    function testHaldanePeakLocation(tc)
    %   The methanogenic growth rate is maximal at S2 = sqrt(KS2*KI2); past
    %   that point more substrate means slower growth, which is the
    %   mechanism behind digester souring.
        S2 = linspace(0.1, 400, 4000)';
        k = am2_kinetics(10*ones(size(S2)), S2, 37, tc.D, tc.P);
        [~, i] = max(k.mu2);
        expected = sqrt(tc.D.am2.KS2 * tc.D.am2.KI2);
        tc.verifyEqual(S2(i), expected, 'RelTol', 0.02);
    end

    function testAM2SteadyStateMatchesIntegration(tc)
    %   Integrate the isolated digester at constant temperature and feed,
    %   and check it lands on the closed-form equilibrium. This validates
    %   the analytical solution used for BMP calibration against the same
    %   mass balances the full model integrates.
        Dd = tc.D;
        Dil  = 1/25;
        S1in = 20;
        S2in = Dd.feed.S2in_mmolL;
        T    = 37;

        ss = am2_steady_state(Dil, S1in, S2in, T, Dd, tc.P);
        tc.assumeFalse(ss.washout, 'Test point must not wash out.');

        x0 = [1; 0.5; 10; 0.5];
        f  = @(t,x) am2_only_rhs(x, Dil, S1in, S2in, T, Dd, tc.P);
        s  = ode15s(f, [0 3000], x0, odeset('RelTol',1e-9,'AbsTol',1e-11));
        xe = s.y(:, end);

        tc.verifyEqual(xe(1), ss.S1, 'RelTol', 1e-3);
        tc.verifyEqual(xe(2), ss.X1, 'RelTol', 1e-3);
        tc.verifyEqual(xe(3), ss.S2, 'RelTol', 1e-3);
        tc.verifyEqual(xe(4), ss.X2, 'RelTol', 1e-3);
    end

    function testWashoutAtHighDilution(tc)
    %   Above the maximum growth rate the populations cannot keep up with
    %   removal and the reactor washes out. The closed form must report
    %   that rather than returning a negative or complex biomass.
        Dil = 5;                 % far above mu2max = 0.74 /d
        ss = am2_steady_state(Dil, 20, 5, 37, tc.D, tc.P);
        tc.verifyTrue(ss.washout);
        tc.verifyEqual(ss.X2, 0);
        tc.verifyEqual(ss.qM, 0);
    end

    function testBMPCalibrationHitsTarget(tc)
    %   k6 appears only in the methane output equation, never in the mass
    %   balances, so the calibration is exactly linear and must land on the
    %   target in one step.
        [~, info] = calibrate_am2_to_bmp(tc.D, tc.P, 25, 40);
        tc.verifyEqual(info.BMP_after, info.BMP_target, 'RelTol', 1e-9);
    end

    function testGasCompositionIsPhysical(tc)
    %   Methane fraction must be a fraction, and pH must land in a range a
    %   digester could actually occupy.
        g = am2_gas(0.6, 50, 55, 0.03, 0.65, 4.3, tc.D, tc.P);
        tc.verifyGreaterThan(g.x_CH4, 0);
        tc.verifyLessThan(g.x_CH4, 1);
        tc.verifyGreaterThan(g.pH, 5);
        tc.verifyLessThan(g.pH, 9);
        tc.verifyGreaterThanOrEqual(g.Q_CH4_Nm3_d, 0);
    end

    function testArrheniusFactorMatchesStaticModel(tc)
    %   The default temperature response must reproduce the exponential the
    %   static model applies to BMP, capped at the reference temperature.
        T = (10:1:45)';
        fT = growth_temp_factor(T, tc.D, tc.P);
        expected = min(exp(tc.P.theta_arrhenius*(T - tc.D.temp.T_ref)), 1);
        tc.verifyEqual(fT, expected, 'AbsTol', 1e-12);
    end

    function testCTMICardinalPoints(tc)
    %   The cardinal temperature model must be 1 at the optimum and vanish
    %   at both cardinal extremes.
        Dd = tc.D;  Dd.temp.model = 'ctmi';
        tc.verifyEqual(growth_temp_factor(Dd.temp.T_opt, Dd, tc.P), 1, ...
                       'AbsTol', 1e-10);
        tc.verifyEqual(growth_temp_factor(Dd.temp.T_min, Dd, tc.P), 0, ...
                       'AbsTol', 1e-10);
        tc.verifyEqual(growth_temp_factor(Dd.temp.T_max, Dd, tc.P), 0, ...
                       'AbsTol', 1e-10);
        tc.verifyEqual(growth_temp_factor(50, Dd, tc.P), 0, 'AbsTol', 1e-12);
    end

    end

% =====================================================================
    methods (Test)  % --- field sub-models ----------------------------

    function testWaterStressCoefficientBounds(tc)
    %   Ks must stay within [0,1] and decrease monotonically as the root
    %   zone dries.
        TAW = 1000*(tc.D.water.theta_fc - tc.D.water.theta_wp)*tc.D.water.Zr_m;
        Dr = linspace(0, 1.5*TAW, 500)';
        Ks = zeros(size(Dr));
        for i = 1:numel(Dr)
            w = soil_water(Dr(i), 5, 1.0, 0, 0, 1, tc.D, tc.P);
            Ks(i) = w.Ks;
        end
        tc.verifyTrue(all(Ks >= 0 & Ks <= 1 + 1e-12));
        tc.verifyTrue(all(diff(Ks) <= 1e-9), 'Ks must be non-increasing in Dr.');
    end

    function testShadingReducesEvapotranspiration(tc)
    %   The shading term must cut actual crop ET by exactly the coefficient
    %   the static model applies, so the two water accountings agree.
        w0 = soil_water(0, 6, 1.0, 0.0, 0, 1, tc.D, tc.P);
        w1 = soil_water(0, 6, 1.0, 0.5, 0, 1, tc.D, tc.P);
        expected = w0.ETc_mm_d * (1 - tc.P.alpha_shade*0.5);
        tc.verifyEqual(w1.ETc_mm_d, expected, 'RelTol', 1e-12);
    end

    function testBiomassIsLinearInRUE(tc)
    %   CALIBRATE_CROP_RUE relies on exact linearity: canopy development
    %   depends on thermal time alone and the harvest term is linear in
    %   biomass, so scaling RUE scales growth by the same factor. If this
    %   ever stops holding, the one-shot calibration silently becomes wrong.
        Dd = tc.D;
        c1 = crop_growth(500, 2000, 200, 174, 22, 1, 0.9, 0, Dd);
        Dd.crop.RUE_g_MJ = 2*tc.D.crop.RUE_g_MJ;
        c2 = crop_growth(500, 2000, 200, 174, 22, 1, 0.9, 0, Dd);
        tc.verifyEqual(c2.growth_kgDM_d, 2*c1.growth_kgDM_d, 'RelTol', 1e-12);
    end

    function testPARSaturationCaps(tc)
    %   Above the saturation point extra light must not produce extra
    %   growth -- the property that makes agrivoltaic shading affordable.
    %   The tolerance is set by the smoothing width of SMOOTH_MIN, not by
    %   floating point: at 2 % of PAR_sat the smoothed cap still creeps by
    %   about 5e-5 in relative terms between 400 and 900 W/m2. That residual
    %   is the deliberate price of a differentiable right-hand side.
        Dd = tc.D;
        cLo = crop_growth(500, 0, 400, 174, 22, 1, 1, 0, Dd);
        cHi = crop_growth(500, 0, 900, 174, 22, 1, 1, 0, Dd);
        tc.verifyEqual(cHi.growth_kgDM_d, cLo.growth_kgDM_d, 'RelTol', 1e-3);

        % A hard cap must be recovered as the smoothing width goes to zero.
        pLo = smooth_min(400, 174, 1e-9);
        pHi = smooth_min(900, 174, 1e-9);
        tc.verifyEqual(pLo, 174, 'AbsTol', 1e-6);
        tc.verifyEqual(pHi, 174, 'AbsTol', 1e-6);
    end

    function testNoGrowthOutOfSeason(tc)
        c = crop_growth(500, 1000, 300, 174, 25, 0, 1, 0, tc.D);
        tc.verifyEqual(c.growth_kgDM_d, 0, 'AbsTol', 1e-12);
    end

    function testThermalTimeResetsBetweenSeasons(tc)
    %   REGRESSION. Thermal time is a state; if it is not reset out of
    %   season, a multi-year run begins its second season on a canopy that
    %   is already past senescence and the crop produces essentially
    %   nothing. The failure is nasty because every RATIO still looks
    %   plausible -- the open-field reference run carries the same
    %   accumulated thermal time, so LER_crop survives while the absolute
    %   yields collapse and the digester quietly starves.
        c = crop_growth(1500, 0, 0, 174, 20, 0, 1, 0, tc.D);
        tc.verifyLessThan(c.dGDD_dt, 0, ...
            'Thermal time must decay when out of season.');
        tau = tc.D.crop.tau_gdd_reset_d;
        tc.verifyEqual(c.dGDD_dt, -1500/tau, 'RelTol', 1e-12);

        % Sixty days out of season -- the shortest gap of the four sites --
        % must clear it essentially completely.
        s = ode45(@(t,g) -g/tau, [0 60], 1500);
        tc.verifyLessThan(s.y(end), 1e-6*1500);
    end

    function testHarvestRemovesTheStandingCrop(tc)
    %   REGRESSION. Removal is first-order in biomass, so a unit-area
    %   harvest pulse leaves exp(-1) = 37 % of the crop standing in the
    %   field. The pulse strength must be large enough for a harvest to be
    %   a harvest.
        tc.verifyGreaterThanOrEqual(tc.D.feed.harvest_strength, 4, ...
            'A harvest must remove at least 98 % of standing biomass.');

        c = tc.getClimate();
        F = build_forcing(c, tc.site, tc.design, tc.P, tc.D, 'apv');
        tq = (0:1/24:365-1/24)';
        tc.verifyEqual(trapz(tq, F.harvest(tq)), ...
                       tc.D.feed.harvest_strength, 'RelTol', 0.02);
    end

    end

% =====================================================================
    methods (Test)  % --- assembled model -----------------------------

    function testRHSIsFiniteAndCorrectlySized(tc)
        c = tc.getClimate();
        F = build_forcing(c, tc.site, tc.design, tc.P, tc.D, 'apv');
        x0 = apvad_initial_state(c, F, tc.D, tc.P);
        ix = apvad_states();
        for t = [0 90 180 270 364]
            dx = apvad_rhs(t, x0, F, tc.D, tc.P);
            tc.verifySize(dx, [ix.n 1]);
            tc.verifyTrue(all(isfinite(dx)), ...
                sprintf('Non-finite derivative at t = %g d.', t));
        end
    end

    function testHeaterIsSizedToTheDesignLoad(tc)
    %   An oversized heater holds the setpoint in every scenario and
    %   destroys the model's ability to say whether the diverted PV
    %   fraction is sufficient. Capacity must track the computed load.
        c = tc.getClimate();
        F = build_forcing(c, tc.site, tc.design, tc.P, tc.D, 'apv');
        tc.verifyEqual(F.Q_max_W, tc.D.thermal.Q_max_factor*F.Q_design_W, ...
                       'RelTol', 1e-12);
        tc.verifyLessThan(F.Q_max_W, 5000, ...
            'A few-cubic-metre farm digester should need well under 5 kW.');
        tc.verifyGreaterThan(F.Q_max_W, 50);
    end

    function testAccumulatorRatesAreNonNegative(tc)
    %   The energy, water and harvest integrators are fed output rates that
    %   are non-negative by construction, so a negative rate anywhere is a
    %   sign error upstream.
    %
    %   This is asserted on the RATES rather than on the integrated series,
    %   and the distinction matters. ODE15S is an implicit multistep method
    %   controlling local error to RelTol; on an accumulator of order 1e4
    %   kWh that permits step-to-step wobble of order 0.1 kWh, so the
    %   integrated series is monotone only to solver tolerance and a strict
    %   monotonicity test on it measures the solver, not the model. The
    %   rates are exact.
        c = tc.getClimate();
        F = build_forcing(c, tc.site, tc.design, tc.P, tc.D, 'apv');
        ix = apvad_states();
        x0 = apvad_initial_state(c, F, tc.D, tc.P);
        opt = odeset('RelTol',1e-5,'AbsTol',1e-7,'MaxStep',1/24, ...
                     'JPattern',ix.JPattern);
        s = ode15s(@(t,x) apvad_rhs(t,x,F,tc.D,tc.P), [0 15], x0, opt);

        tq = linspace(0, 15, 600);
        Xq = deval(s, tq);
        for i = 1:numel(tq)
            dx = apvad_rhs(tq(i), Xq(:,i), F, tc.D, tc.P);
            for j = ix.accumulators
                tc.verifyGreaterThanOrEqual(dx(j), -1e-12, ...
                    sprintf('Rate of %s negative at t = %.3f d.', ...
                            ix.names{j}, tq(i)));
            end
        end

        % Over the whole window the integrals must still grow.
        for j = ix.accumulators
            tc.verifyGreaterThanOrEqual(s.y(j,end) - s.y(j,1), -1e-9, ...
                sprintf('Accumulator %s ended below its start.', ix.names{j}));
        end
    end

    function testReferenceModeRemovesGroundShading(tc)
    %   The open-field reference supplies every LER denominator, so it must
    %   genuinely carry no inter-row ground shading.
        c = tc.getClimate();
        Fr = build_forcing(c, tc.site, tc.design, tc.P, tc.D, 'reference');
        tc.verifyEqual(max(Fr.F_shad(Fr.t_day)), 0, 'AbsTol', 1e-12);
        Fa = build_forcing(c, tc.site, tc.design, tc.P, tc.D, 'apv');
        tc.verifyGreaterThan(max(Fa.F_shad(Fa.t_day)), 0.1);
    end

    function testForcingIsPeriodic(tc)
    %   Multi-year runs repeat the same typical year, so the forcing must
    %   wrap cleanly at the period boundary.
        c = tc.getClimate();
        F = build_forcing(c, tc.site, tc.design, tc.P, tc.D, 'apv');
        tc.verifyEqual(F.T_amb(mod(370, F.t_period)), F.T_amb(370 - 365), ...
                       'AbsTol', 1e-9);
        tc.verifyEqual(F.t_period, 365, 'AbsTol', 1e-6);
    end

    end
end

% =====================================================================
% Local helpers used by the tests
% =====================================================================
function dT = pv_thermal_rhs(T, G, Ta, WS, D, P)
pv = pv_module(G, Ta, WS, T, 1, D, P);
dT = pv.dTcell_dt;
end

function dx = am2_only_rhs(x, Dil, S1in, S2in, T, D, P)
%AM2_ONLY_RHS The four AM2 balances in isolation, at fixed temperature.
%   Z and C are held at their influent values, which is their exact steady
%   state when the dilution rate is constant.
a  = D.am2;
S1 = max(x(1),0); X1 = max(x(2),0); S2 = max(x(3),0); X2 = max(x(4),0);
k  = am2_kinetics(S1, S2, T, D, P);
dx = [ Dil*(S1in - S1) - a.k1*k.mu1*X1
       (k.mu1 - a.alpha*Dil - a.kd1)*X1
       Dil*(S2in - S2) + a.k2*k.mu1*X1 - a.k3*k.mu2*X2
       (k.mu2 - a.alpha*Dil - a.kd2)*X2 ];
end
