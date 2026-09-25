function ix = apvad_states()
%APVAD_STATES Index map, names and units of the coupled state vector.
%
%   IX = APVAD_STATES() returns a struct whose fields are the indices of
%   each state in the vector integrated by ODE15S, plus bookkeeping fields
%   IX.n (number of states), IX.names, IX.units and IX.labels.
%
%   Using named indices everywhere (x(ix.S1) rather than x(1)) means the
%   state vector can be extended without hunting for magic numbers.
%
%   State vector (all rates are per DAY; t = 0 is 1 January 00:00)
%   ------------------------------------------------------------------
%   Anaerobic digestion -- AM2 (Bernard et al. 2001)
%    1  S1      g COD/L    soluble/hydrolysable organic substrate
%    2  X1      g VSS/L    acidogenic biomass
%    3  S2      mmol/L     volatile fatty acids
%    4  X2      g VSS/L    methanogenic biomass
%    5  Z       mmol/L     total alkalinity
%    6  C       mmol/L     total inorganic carbon
%
%   Thermal
%    7  T_dig   degC       digester bulk temperature
%    8  T_cell  degC       PV module temperature
%
%   Field
%    9  Dr      mm         root-zone water depletion (FAO-56)
%   10  GDD     degC d     accumulated growing degree days
%   11  B_dm    kg DM/ha   standing above-ground crop dry biomass
%   12  M_res   kg VS      feedstock stored and available to the digester
%
%   Accumulators (pure integrators of output rates)
%   13  E_pv        kWh/ha  cumulative PV generation
%   14  E_ch4       kWh/ha  cumulative methane energy produced
%   15  E_heat      kWh/ha  cumulative heat delivered to the digester
%   16  E_pv_sold   kWh/ha  cumulative PV electricity exported
%   17  E_ch4_heat  kWh/ha  cumulative methane burned in the backup boiler
%   18  W_irr       mm      cumulative irrigation applied
%   19  M_harv      kg DM/ha cumulative harvested above-ground dry matter
%
%   Carrying the accumulators inside the ODE means annual totals come
%   straight from the solver, under its own error control, instead of being
%   re-integrated from a resampled output grid. They are linear and add no
%   stiffness.
%
%   The result is cached in a PERSISTENT variable. APVAD_RHS calls this on
%   every evaluation, and the solver evaluates the right-hand side hundreds
%   of thousands of times per simulated year; rebuilding three cell arrays
%   and nineteen struct fields each time dominated the runtime before the
%   cache was added. Call CLEAR APVAD_STATES if the state vector is ever
%   edited during a session.
%
%   See also APVAD_RHS, APVAD_INITIAL_STATE.

persistent cached
if ~isempty(cached)
    ix = cached;
    return
end

names  = {'S1','X1','S2','X2','Z','C', ...
          'T_dig','T_cell','Dr','GDD','B_dm','M_res', ...
          'E_pv','E_ch4','E_heat','E_pv_sold','E_ch4_heat','W_irr','M_harv'};

units  = {'g COD/L','g VSS/L','mmol/L','g VSS/L','mmol/L','mmol/L', ...
          'degC','degC','mm','degC d','kg DM/ha','kg VS', ...
          'kWh/ha','kWh/ha','kWh/ha','kWh/ha','kWh/ha','mm','kg DM/ha'};

labels = {'Organic substrate S_1','Acidogenic biomass X_1', ...
          'Volatile fatty acids S_2','Methanogenic biomass X_2', ...
          'Total alkalinity Z','Inorganic carbon C', ...
          'Digester temperature','PV module temperature', ...
          'Root-zone depletion','Growing degree days', ...
          'Standing crop biomass','Stored feedstock', ...
          'Cumulative PV generation','Cumulative methane energy', ...
          'Cumulative digester heat','Cumulative PV exported', ...
          'Cumulative methane to boiler','Cumulative irrigation', ...
          'Cumulative harvested dry matter'};

for k = 1:numel(names)
    ix.(names{k}) = k;
end
ix.n      = numel(names);
ix.names  = names;
ix.units  = units;
ix.labels = labels;

% Indices of the pure accumulators, used by the postprocessor.
ix.accumulators = [ix.E_pv, ix.E_ch4, ix.E_heat, ix.E_pv_sold, ...
                   ix.E_ch4_heat, ix.W_irr, ix.M_harv];

% Jacobian sparsity pattern. Nothing in the model reads an accumulator --
% they are write-only integrators of output rates -- so their COLUMNS are
% structurally zero. Handing this to ODE15S lets NUMJAC group and skip
% those columns, cutting the finite-difference cost of every Jacobian by
% about a third. The remaining block is genuinely dense: the digester
% temperature reaches the kinetics, the kinetics reach the gas phase, the
% gas phase reaches the boiler and back to the temperature.
ix.JPattern = sparse(ones(ix.n, ix.n));
ix.JPattern(:, ix.accumulators) = 0;

cached = ix;
end
