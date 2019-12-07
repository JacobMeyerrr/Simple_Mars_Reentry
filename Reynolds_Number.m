function Re = Reynolds_Number(v, D, nu)
%% Reynolds_Number
%  Computes the reynolds number for flow over the spacecraft.
%
% USAGE:
%       Reynolds_Number(v,D,nu)
%
% INPUTS:
%       v: Magnitude of velocity for the spacecraft                 (km/s)
%       D: Diameter of spacecraft leading edge                      (m)
%      nu: Kinematic viscosity for current atmospheric conditions   (m^2/s)
%
% OUTPUTS:
%      Re: Reynolds Number for the flow over the spacecraft      (unitless)
%
%% Function Main

% conver velocity units to (m/s)
v = v*1000;

%calculate Reynolds Number

Re = v*D./nu;

