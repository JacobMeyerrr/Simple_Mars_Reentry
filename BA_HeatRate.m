function qAvg = BA_HeatRate(rho, v ,Cf)
%% BA_HeatRate
%  Calculates the Body-Average heating rate for the spacecraft.
%
% USAGE:
%       BA_HeatRate(rho, v, Cf)
%
% INPUTS:
%       rho:  The current atmospheric density               (kg/m^3)
%         v:  The spacecrafts current velocity magnitude    (km/s)
%        Cf:  The body-averaged skin coefficient            (unitless)
%
% OUTPUT:
%      qAvg:  The body-averaged heating rate for the spacecraft (kg/s^3)
%
%% Function Main
% Convert km/s to m/s
v = v*1000;

% Calculate Body-Averaged Heating Rate
qAvg = 1/4*rho.*v.^3.*Cf;

