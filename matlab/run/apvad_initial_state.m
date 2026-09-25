function x0 = apvad_initial_state(clim, F, D, P)
%APVAD_INITIAL_STATE Physically plausible initial condition for the model.
%
%   X0 = APVAD_INITIAL_STATE(CLIM, F, D, P) returns the initial state
%   vector at t = 0 (1 January 00:00).
%
%   The digester is initialised as an established, seeded, buffered reactor
%   rather than as a cold empty vessel: a start-up transient from zero
%   biomass takes months and would contaminate the first simulated year.
%   RUN_DYNAMIC_SIM additionally discards a full spin-up year by default,
%   so the reported results depend on this initial guess only through
%   whether it lands in the basin of the healthy operating point rather
%   than the acidified one -- which is itself a result worth reporting, and
%   is what TESTWASHOUTATHIGHLOADING probes.
%
%   The field starts the calendar year at field capacity (Dr = 0), out of
%   season, with no standing biomass, and with a feedstock store holding
%   the residues of the previous season.
%
%   See also RUN_DYNAMIC_SIM, APVAD_STATES.

ix = apvad_states();
x0 = zeros(ix.n, 1);

% ---- Digester: seeded and buffered -----------------------------------
x0(ix.S1) = 1.0;      % g COD/L   low residual substrate
x0(ix.X1) = 0.5;      % g VSS/L   established acidogens
x0(ix.S2) = 10.0;     % mmol/L    modest VFA pool, well below inhibition
x0(ix.X2) = 0.5;      % g VSS/L   established methanogens
x0(ix.Z)  = D.feed.Zin_mmolL;
x0(ix.C)  = D.feed.Cin_mmolL;

% ---- Thermal ---------------------------------------------------------
T_amb_mean = mean(clim.T2m);
x0(ix.T_dig)  = min(max(T_amb_mean + P.T_amb_boost_degC, ...
                        P.T_dig_min_degC), P.T_ref_degC);
x0(ix.T_cell) = clim.T2m(1);

% ---- Field -----------------------------------------------------------
x0(ix.Dr)    = 0;     % mm        at field capacity on 1 January
x0(ix.GDD)   = 0;
x0(ix.B_dm)  = 0;

% Carry-over feedstock, started at the store's own equilibrium level:
%     M* = buffer_days * (daily feed demand)
% so that dM/dt = supply - feed - k*M is satisfied from the first day. The
% store's time constant is 1/k = 1000 days, so starting anywhere else means
% the simulation spends decades drifting towards this value and reports
% every intervening year as if it were a steady one. See the note on
% D.feed.buffer_days in APVAD_DYN_PARAMS.
x0(ix.M_res) = D.feed.buffer_days * F.Q_in_m3_d * F.c_VS_in_gL;

% ---- Accumulators start at zero --------------------------------------
x0(ix.accumulators) = 0;
end
