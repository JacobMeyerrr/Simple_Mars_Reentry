% Computes the gravitational force of Mars on a spacecraft in orbit or the
% atmosphere
%
% Usage:
%   [g] = gravityforce(r,m)
%
% Inputs:
%   r = radius vector around mars (km)
%   m = mass of your spacecraft (kg)
%
% Outputs:
%   g = the gravity vector [gx gy gz] of current state (m/s/s)
%

function [g] = gravityforce(r,m)

% define Mars gravitational constant
mu = 4.282837E13; % (km^3/s^2)

% convert r to meters instead of km
r = r*1000;

% obtain magnitude of r 
R = sqrt(sum(r.^2));

% compute final result (gravity)
g = (-m*mu/R^3).*r

end


