function f = feed_scheduler(M_res, residue_VS_in, F, D)
%FEED_SCHEDULER Draw feedstock from storage and set the influent composition.
%
%   F_OUT = FEED_SCHEDULER(M_RES, RESIDUE_VS_IN, F, D) returns
%
%     F_OUT.VS_feed_kg_d   volatile solids actually fed (kg VS/d)
%     F_OUT.c_VS_in_gL     influent VS concentration (g VS/L)
%     F_OUT.S1in           influent organic substrate (g COD/L)
%     F_OUT.S2in           influent VFA (mmol/L)
%     F_OUT.Zin, F_OUT.Cin influent alkalinity and inorganic carbon (mmol/L)
%     F_OUT.dM_res_dt      storage derivative (kg VS/d)
%     F_OUT.OLR_kgVS_m3d   instantaneous organic loading rate
%     F_OUT.avail          storage availability factor (0..1)
%
%   The digester is fed at a constant volumetric rate Q_in = V/HRT. What
%   varies is the influent concentration, because feedstock arrives on the
%   farm's calendar and not the digester's: cattle manure trickles in all
%   year, while the entire tomato residue crop lands within a few days of
%   harvest. The store buffers that mismatch, and its dynamics are:
%
%       dM_res/dt = residue_in + manure - VS_feed - storage_loss*M_res
%
%   The nominal feed demand Q_in*c_VS_in is set in BUILD_FORCING so that
%   the annual average loading rate equals the static model's
%   OLR = VS_annual/(365*V_dig). When the store runs dry the availability
%   factor closes the draw smoothly and the digester is under-fed -- which
%   is the honest outcome, and one the static model cannot express: dividing
%   an annual VS total by 365 days implicitly assumes a feedstock store of
%   unlimited size that is always exactly full.
%
%   See also BUILD_FORCING, AM2_KINETICS, SMOOTH_STEP.

Q_in = F.Q_in_m3_d;

% ---- Nominal demand and store availability ---------------------------
% The availability switch closes over half a day's feed demand rather than
% over a fixed kilogram, so the transition is smooth relative to the rate
% at which the store is actually drawn down.
VS_demand = Q_in * F.c_VS_in_gL;                       % kg VS/d
sw        = max([D.feed.min_store_kgVS, 0.5*VS_demand, eps]);
avail     = smooth_step(M_res - D.feed.min_store_kgVS, sw);
VS_feed   = VS_demand * avail;

% ---- Influent composition --------------------------------------------
c_VS = VS_feed / max(Q_in, eps);                       % kg/m3 = g/L
S1in = D.feed.f_S1 * D.feed.COD_per_VS * c_VS;         % g COD/L

% ---- Storage balance --------------------------------------------------
dM_res_dt = residue_VS_in + F.manure_VS_kg_d - VS_feed ...
            - D.feed.store_loss_1_d * max(M_res, 0);

f.VS_feed_kg_d = VS_feed;
f.c_VS_in_gL   = c_VS;
f.S1in         = S1in;
f.S2in         = D.feed.S2in_mmolL;
f.Zin          = D.feed.Zin_mmolL;
f.Cin          = D.feed.Cin_mmolL;
f.dM_res_dt    = dM_res_dt;
f.OLR_kgVS_m3d = VS_feed / max(F.V_dig_m3, eps);
f.avail        = avail;
end
