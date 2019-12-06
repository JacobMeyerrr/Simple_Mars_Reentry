function StateDeriv = CraftOrbit(time, State, Mbody, Mcraft, G)
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
%      Mbody: Mass of the parent body. Units: [kg]
%     Mcraft: Mass of the spacecraft.  Units: [kg]
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

Rhat = [x;y;z]/R;

%Velocity
xDot = State(4);            %(km/s)
yDot = State(5);            %(km/s)
zDot = State(6);            %(km/s)

V = sqrt(xDot^2+yDot^2+zDot^2);

Vhat = [xDot; yDot; zDot]/V;

%Forces
Fg = -G*Mbody*Mcraft/(R*1000)^2;   %(N)

Fd = 0;

%Acceleration
ax = (Fg*Rhat(1)/1000) ...
             + (Fd*-Vhat(1));       %kg*km/s^2
            
ay = (Fg*Rhat(2)/1000) ...
             + (Fd*-Vhat(2));       %kg*km/s^2
            
az = (Fg*Rhat(3)/1000) ...
             + (Fd*-Vhat(3));       %kg*km/s^2
            

StateDeriv = [xDot; yDot; zDot; ax; ay; az];

