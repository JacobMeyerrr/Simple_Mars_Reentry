
%Define constants
Radius_Mars = 3390;                     %(km)
Mass_Mars   = 6.39e23;                  %(kg)
Mass_craft  = 10400;                    %(kg)
G           = (6.67408e-11)/1000^3;     %(km^3/(kg*s^2))

%Define initial conditions
Altitude = Radius_Mars + 200;           %(km)
V0 = sqrt(G*(Mass_Mars)/Altitude);                              %(km/s)

%Define initial state vector and time vector
PosState = [Altitude;0;0];              %[x;y;z]   (km;km;km)
VelState = [0;V0;0];                    %[Vx;Vy;Vz](km/s;km/s;km/s;)

state0 = [PosState; VelState];

time = linspace(1,6600,9e4);

%Calculating orbit
orbit = ODENumIntRK4(@CraftOrbit,time,state0,Mass_Mars,Mass_craft,G);


figure(1), hold on, grid on
plot3(orbit(:,1),orbit(:,2),orbit(:,3),'b-')
pbaspect([1 1 1])


