function M = Mach_Number(v, s)
%% Mach_Number
%  Computes the mach number of the spacecraft for its current velocity and
%  altitude in the Martian atmosphere
%
% USAGE:
%       Mach_Number(v, s)
%
% INPUTS:
%       v: Magnitude of the spacecrafts current velocity (km/s)
%       s: speed of sound in the current for the current atmospheric
%          conditions (m/s)
%
% OUTPUTS:
%       M: Mach number (unitless) 
%
%% Function Main

M = v/s;