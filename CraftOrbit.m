function [StateDeriv] = CraftOrbit(time, State, Mcraft, Drag0_1)
%% Craft Orbit
%  Calculate the orbital parameters of the spacecraft in orbit of the
%  parent body
%  
%  Right-hand side for the dynamic system of a spacecraft in orbit
%
% USAGE:
%       StateDeriv = CraftOrbit(time, State, Mbody, Mcraft)
%
% INPUTS:
%       time: Current time. (Not used in this calculation)
%      State: State vector at current time.
%             [x; y; z; Vx; Vy; Vz] - Units: [m; m; m; m/s; m/s; m/s]
%
%     Mcraft: Mass of the spacecraft.  Units: [kg]
%    Drag0_1: Boolean input to turn drag on/off (true = on)
%
% OUTPUTS:
%        StateDeriv = [xDot; yDot; zDot; aX; aY; aZ] 
%                     Units: [m/s; m/s; m/s; m/s^2; m/s^2; m/s^2]
%

%% Function Main

%Postition information
x = State(1);               %(km)
y = State(2);               %(km)
z = State(3);               %(km)

R = sqrt(x^2+y^2+z^2);      %(km)

%Velocity
xDot = State(4);            %(km/s)
yDot = State(5);            %(km/s)
zDot = State(6);            %(km/s)

V = sqrt(xDot^2+yDot^2+zDot^2);

Vhat = [xDot; yDot; zDot]/V;

%Calculate density
rho = AtmDensityMars(R-3390);

% Calculate Gravitational Force
Fg = gravityforce([x;y;z],Mcraft);  %(N)

% Calculate Drag Force
if(Drag0_1)
    Fd = DragForce(rho,10.4,20,V);
else
    Fd = 0;
end

% Calculate Acceleration components
ax = (Fg(1)/(Mcraft*1000)) ...
             + (Fd*-Vhat(1)/(Mcraft*1000));       %km/s^2
            
ay = (Fg(2)/(Mcraft*1000)) ...
             + (Fd*-Vhat(2)/(Mcraft*1000));       %km/s^2
            
az = (Fg(3)/(Mcraft*1000)) ...
             + (Fd*-Vhat(3)/(Mcraft*1000));       %km/s^2
            

StateDeriv = [xDot; yDot; zDot; ax; ay; az];

