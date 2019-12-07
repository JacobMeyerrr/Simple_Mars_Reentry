
%Define Planet constants
Radius_Mars = 3390;                     %(km)
Mass_Mars   = 6.39e23;                  %(kg)

%Define Spacecraft Constants
Mass_craft = 10400;                     %(kg)
X_Area     = 19.635;                    %(m) 

G           = (6.67408e-11)/1000^3;     %(km^3/(kg*s^2))

%Define initial conditions
Altitude = Radius_Mars + 200;           %(km)
V0 = sqrt(G*(Mass_Mars)/Altitude);                              %(km/s)

%Define initial state vector and time vector
PosState = [Altitude;0;0];              %[x;y;z]   (km;km;km)
VelState = [0;V0;0];                    %[Vx;Vy;Vz](km/s;km/s;km/s;)

state0 = [PosState; VelState];

time = linspace(1,30000,9e5);          % seconds

%Calculating orbit
orbit = ODENumIntRK4(@CraftOrbit,time,state0,Mass_Mars,Mass_craft,G);


figure(1), hold on, grid on
plot(orbit(:,1),orbit(:,2),'b-')
plot(x,y,'r-')
plot(x,-y,'r-')
pbaspect([1 1 1])


