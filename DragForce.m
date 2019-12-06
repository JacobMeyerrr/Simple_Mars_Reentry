function drag = DragForce(rho,B,m,A,v)

%Drag Force 
%Calculates the drag force of the spacecraft depending on the current
%atmospheric density, ballistic coefficient, mass, cross-sectional area,
%and velocity.
%
%Usage:
%   DragForce(rho,B,m,A,v)
%
%Inputs:
%   rho: atmospheric density    [kg/m^3]
%   B: ballistic coefficient    [kg/m^2]
%   m: mass                     [kg]
%   A: cross-sectional area     [m^2]
%   v: velocity                 [km/s]
%
%Outputs:
%   drag: drag force            [N]         

v=v*1000;

Cd=m/(B*A);

drag=1/2*rho*Cd*A*v.^2;

