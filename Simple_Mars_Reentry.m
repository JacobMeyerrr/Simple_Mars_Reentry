
%Define Planet constants
Radius_Mars = 3390;                     %[km]
Mass_Mars   = 6.39e23;                  %[kg]

%Define Spacecraft Constants
Mass_craft = 10400;                     %[kg]
X_Area     = 19.635;                    %[m]

G          = (6.67408e-11)/1000^3;     %[km^3/(kg*s^2)]

%Define initial conditions
Altitude = Radius_Mars + 200;           %[km]
V0_real = sqrt(G*(Mass_Mars+Mass_craft)/Altitude)+0.0073062;     %[km/s]
V0 = V0_real*0.5;

%Define initial state vector and time vector

PosState = [Altitude;0;0];              %[x;y;z]   [km;km;km]
VelState = [0;V0;0];                    %[Vx;Vy;Vz][km/s;km/s;km/s]

state0 = [PosState; VelState];


T = sqrt(4*pi^2*(norm(PosState))^3/(G*Mass_Mars));  %calculating the orbital period
NumOrbs = 50;
time = linspace(0,T*NumOrbs,1e5);          % [s]

%Calculating orbit using numerical integration
orbit = ODENumIntRK4(@CraftOrbit,time,state0,Mass_craft);

figure(1), hold on, grid on
pbaspect([1 1 1])
plot(x,y,'r-')
plot(x,-y,'r-')
plot(orbit(:,1),orbit(:,2),'b-')

%Calculate Body-Averaged HeatingRate
qAvg = HeatingRate(orbit);

R = vecnorm(orbit(:,1:3),2,2)-Radius_Mars;

figure(2), grid on
plot(R,qAvg,'c-')

% 7. Calculate Max Heating Rate, and Alt where it occurs
maxQ = max(qAvg)
alt_MaxQ = R(find(qAvg == maxQ))

% 8, Record Final Velocity, Flight Path Angle (Gamma), and time elapsed
vfinal = norm(orbit(end,4:6))
gamma_final = flightpathangle(orbit(end,1:3),orbit(end,4:6))
elapsed_time = sum(time(1,length(orbit)))




