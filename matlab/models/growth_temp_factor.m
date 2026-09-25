function fT = growth_temp_factor(T_degC, D, P)
%GROWTH_TEMP_FACTOR Dimensionless temperature response of microbial growth.
%
%   FT = GROWTH_TEMP_FACTOR(T_DEGC, D, P) returns the factor multiplying
%   the maximum specific growth rates mu1max and mu2max. Two forms are
%   available, selected by D.temp.model:
%
%   'arrhenius' (default)
%       fT = exp(theta * (T - T_ref)),  capped at 1 for T > T_ref
%
%       This is the exponential form the static Python model applies to
%       BMP (PARAMS.theta_arrhenius = 0.069 1/degC, Pilarski 2025). Keeping
%       it as the default means the dynamic model's temperature sensitivity
%       is identical to the static one, so any difference in annual methane
%       output is attributable to the dynamics rather than to a changed
%       temperature law. Its weakness is that it cannot represent thermal
%       inhibition above the optimum.
%
%   'ctmi'
%       Cardinal temperature model with inflection (Rosso et al. 1993):
%
%           fT = (T-Tmax)(T-Tmin)^2 /
%                [ (Topt-Tmin) * ( (Topt-Tmin)(T-Topt)
%                                - (Topt-Tmax)(Topt+Tmin-2T) ) ]
%
%       zero outside [Tmin, Tmax], equal to 1 at Topt. Use this when the
%       digester can exceed the mesophilic optimum -- at Ouagadougou the
%       static model's own effective digester temperature already sits at
%       33 degC, so summer excursions past 40 degC are plausible and the
%       Arrhenius form would reward them instead of penalising them.
%
%   Reference
%     Rosso, L., Lobry, J.R., Flandrois, J.P. (1993). J. Theor. Biol.,
%       162(4), 447-463.
%
%   See also AM2_KINETICS, APVAD_DYN_PARAMS.

% Fast path first: this function is on the innermost loop of the solver, so
% the default case avoids both LOWER and the SWITCH dispatch.
if strcmp(D.temp.model, 'arrhenius')
    fT = exp(P.theta_arrhenius * (T_degC - D.temp.T_ref));
    fT = min(max(fT, 0.0), 1.0);
    return
end

switch lower(D.temp.model)
    case 'arrhenius'
        fT = exp(P.theta_arrhenius * (T_degC - D.temp.T_ref));
        fT = min(max(fT, 0.0), 1.0);

    case 'ctmi'
        Tmin = D.temp.T_min;
        Topt = D.temp.T_opt;
        Tmax = D.temp.T_max;

        num = (T_degC - Tmax) .* (T_degC - Tmin).^2;
        den = (Topt - Tmin) .* ( (Topt - Tmin).*(T_degC - Topt) ...
                               - (Topt - Tmax).*(Topt + Tmin - 2*T_degC) );

        fT = zeros(size(T_degC));
        ok = (T_degC > Tmin) & (T_degC < Tmax) & (abs(den) > eps);
        fT(ok) = num(ok) ./ den(ok);
        fT = min(max(fT, 0), 1);

    otherwise
        error('growth_temp_factor:unknownModel', ...
              'D.temp.model must be ''arrhenius'' or ''ctmi'' (got ''%s'').', ...
              D.temp.model);
end
end
