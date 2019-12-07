function Cf = BA_SkinFric(Re, M)
%% BA_SkinFric
%  Computes the the Body-Average skin friction coefficient.
%
% USAGE:
%       BA_SkinFric(Re, M)
%
% INPUTS:
%       Re: Reynolds number for current atmospheric entry regime (unitless)
%        M: Mach Number for the spacecraft                       (unitless)
%
% OUTPUTS:
%       Cf: Body-Average skin friction coefficient               (unitless)
%
%% Function Main

Cf = (0.65+0.339*((2/pi)*atan(10-M)+1))./sqrt(Re);