
%% Mars Orbital Reentry Simulation
% Simulate the orbital path around Mars for the Nasa Orion capsual.
% Implement deorbit mavuever and simulate entry into the Martian atmosphere
% and terminate the simulation at 5km.
%
% Calculate the Body-Averaged heat rate during reentry
% Calculate the final velocity, flightpath angle, and elapsed time.

%% Defining Constants

% Planet constants
Radius_Mars = 3390;                     %[km]
Mass_Mars   = 6.39e23;                  %[kg]

% Spacecraft Constants
Mass_craft = 10400;                     %[kg]
Diameter   = 5.05;                      %[m]

% Gravitational Constant
G          = (6.67408e-11)/1000^3;      %[km^3/(kg*s^2)]

% Initial conditions
Altitude = Radius_Mars + 200;           %[km]
V0 = sqrt(G*(Mass_Mars+Mass_craft) ... 
              /Altitude)+0.0073062;     %[km/s]

%Define initial state vector and time vector

PosState = [Altitude;0;0];              %[x;y;z]   [km;km;km]
VelState = [0;V0;0];                    %[Vx;Vy;Vz][km/s;km/s;km/s]

state0 = [PosState; VelState];


T = sqrt(4*pi^2*(norm(PosState))^3 ...
                    /(G*Mass_Mars));    %calculating the orbital period
NumOrbs = 3;                            %approximate number of orbits
time1 = linspace(0,T*NumOrbs,1e6);      % [s]

%Pre allocate arrays sizes
maxQ         = zeros(1,10);
alt_MaxQ     = zeros(1,10);
vfinal       = zeros(1,10);
gamma_final  = zeros(1,10);
elapsed_time = zeros(1,10);


%DeltaV's to be tested
DV = linspace(0.9956,0.5,10);


for i = 1:10

% Calculating orbit using numerical integration
orbit = ODENumIntRK4(@CraftOrbit,time1,state0,Mass_craft);

% Initialize deorbit burn
v = orbit(end,4:6);             % current velocity state
u = norm(v);                    % magnitude of velocity
vHat = v/u;                     % velocity unit vector

deltaV = DV(i);                 % DeltaV to induce reentry

Vnew = u*deltaV*vHat;           % New Velocity vector after DeltaV impulse

orbit(end,4:6) = Vnew;          % Orbit(end,:) is our new state vector

time2 = linspace(0,1e6,10e6);   % Define new time vector

% Simulate spacecraft reentry into Mars atmosphere
reentry = ODENumIntRK4(@CraftOrbit,time2,orbit(end,:),Mass_craft);

orbit_Tot = [orbit;reentry];

x = [-3390:1e9:3390];
y = sqrt(3390^2 - x.^2);

figure(i), hold on, grid on
pbaspect([1 1 1])
plot(x,y,'r-')
plot(x,-y,'r-')
plot(orbit_Tot(:,1),orbit_Tot(:,2))

%Calculate Body-Averaged HeatingRate
qAvg = HeatingRate(orbit_Tot);

R = vecnorm(orbit_Tot(:,1:3),2,2)-Radius_Mars;

figure(11), grid on, hold on
plot(R,qAvg)

% Calculate Max Heating Rate, and Alt where it occurs
maxQ(i) = max(qAvg);
alt_MaxQ(i) = R(find(qAvg == maxQ(i)));

% Record Final Velocity, Flight Path Angle (Gamma), and time elapsed
vfinal(i) = norm(orbit_Tot(end,4:6));
gamma_final(i) = flightpathangle(orbit_Tot(end,1:3),orbit_Tot(end,4:6));
elapsed_time(i) = T*NumOrbs + time2(length(orbit_Tot(:,1)))*10;

end

DV = DV*V0 - V0

figure(11) , hold on
leg = legend([num2str(DV(1))], [num2str(DV(2))],[num2str(DV(3))],...
        [num2str(DV(4))],[num2str(DV(5))],[num2str(DV(6))],...
        [num2str(DV(7))],[num2str(DV(8))],[num2str(DV(9))],...
        [num2str(DV(10))]);
title(leg,"deltaV (km/s)") 

format short
t = table(DV(:),elapsed_time(:),vfinal(:),gamma_final(:),maxQ(:), alt_MaxQ(:));
t.Properties.VariableNames = {'DeltaV','Duration','Final_Velocity',...
    'Flight_Path_Angle','Max_Heating_Rate', 'Altitude_MHR'};
t

% figure(12)
% uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
%     'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);







