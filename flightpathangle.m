% Computes the flight path angle of a spacecraft in Mars orbit or
% atmosphere
%
% Usage:
%   [gamma] = flightpathangle(
%
% Inputs:
%   r = radius vector around mars (km)
%   v = velocity vector around mars (m/s/s)
%
% Outputs:
%   gamma = the flight path angle between the horizontal and velocity (rad)
%

function [gamma] = flightpathangle(r,v)

% define the local vertical direction (up) as r_hat
e_up = r/sqrt(sum(r.^2));

% velocity unit vector
e_vel = v/sqrt(sum(r.^2));

% dot product between vertical and velocity unit vectors
d = dot(e_up,e_vel);

% angle between the local vertical and the velocity direction
alpha = acos(d);

% flight path angle
gamma = pi/2 - alpha;

end