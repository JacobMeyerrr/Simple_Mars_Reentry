
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

% Gravitational Constant
G          = (6.67408e-11)/1000^3;      %[km^3/(kg*s^2)]

% Initial conditions
Altitude = Radius_Mars + 200;           %[km]
V0 = sqrt(G*(Mass_Mars + Mass_craft) ... 
              /Altitude)+0.0073062;     %[km/s]

%Define initial state vector and time vector
vHat     = [0;1/sqrt(2);1/sqrt(2)];
PosState = [Altitude;0;0];              %[x;y;z]   [km;km;km]
VelState = V0*vHat;                     %[Vx;Vy;Vz][km/s;km/s;km/s]

state0 = [PosState; VelState];


T = sqrt(4*pi^2*(norm(PosState))^3 ...
                    /(G*Mass_Mars));    %calculating the orbital period
                
NumOrbs = 3;                            %approximate number of orbits
time1 = linspace(0,T*NumOrbs,5e5);      % [s]

%% Calculating the deltaV required to induce reentry
% Define a desired periapsis that will induce a reentry.
% run simulation with varying deltaV's and a zero drag force model
% measuring the periapsis.
% exit simulation when the desired periapsis is achieved
%

dv = linspace(1,0.5,250);     %percent of circular orbital velocity (V0)
entry = 50;                   %Desired periapsis [km]
enters = 0;

for i = 1:length(dv)

    V = V0*dv(i);

    vHat     = [0;1/sqrt(2);1/sqrt(2)];
    PosState = [Altitude;0;0];              %[x;y;z]   [km;km;km]
    VelState = V*vHat;                      %[Vx;Vy;Vz][km/s;km/s;km/s]

    state0 = [PosState; VelState];          %initial state vector

    orbit = ODENumIntRK4(@CraftOrbit,time1,state0,Mass_craft,false);

    R = vecnorm(orbit(:,1:3),2,2)-Radius_Mars;

    periapsis = min(R);

    if(periapsis <= entry)
        initialDV = dv(i);
        break;
    end
    
    if(periapsis <= entry*1.2 && enters == 0)
        dv = linspace(dv(i),dv(i)/2,1e4);
        i=1;
        enters = 1;
    end
end


%% Simulate reentry into Martian atmosphere

% Initial conditions
Altitude = Radius_Mars + 200;           %[km]
V0 = sqrt(G*(Mass_Mars+Mass_craft) ... 
              /Altitude)+0.0073062;     %[km/s]

%Define initial state vector and time vector
vHat     = [0;1/sqrt(2);1/sqrt(2)];
PosState = [Altitude;0;0];              %[x;y;z]   [km;km;km]
VelState = V0*vHat;                     %[Vx;Vy;Vz][km/s;km/s;km/s]

state0 = [PosState; VelState];

%Pre allocate arrays sizes
maxQ         = zeros(1,10);
alt_MaxQ     = zeros(1,10);
vfinal       = zeros(1,10);
gamma_final  = zeros(1,10);
elapsed_time = zeros(1,10);


%DeltaV's to be tested
DV = linspace(initialDV, initialDV/2,10);


for i = 1:10

% Calculating orbit using numerical integration
orbit = ODENumIntRK4(@CraftOrbit,time1,state0,Mass_craft,true);

% Initialize deorbit burn
v = orbit(end,4:6);             % current velocity state
u = norm(v);                    % magnitude of velocity
vHat = v/u;                     % velocity unit vector

deltaV = DV(i);                 % DeltaV to induce reentry

Vnew = u*deltaV*vHat;           % New Velocity vector after DeltaV impulse

orbit(end,4:6) = Vnew;          % Orbit(end,:) is our new state vector

time2 = linspace(0,1e6,1e6);   % Define new time vector

% Simulate spacecraft reentry into Mars atmosphere
reentry = ODENumIntRK4(@CraftOrbit,time2,orbit(end,:),Mass_craft,true);

orbit_Tot = [orbit;reentry];

figure(i), hold on, grid on
[x,y,z] = sphere(50);
surf(x*Radius_Mars,y*Radius_Mars,z*Radius_Mars,'FaceColor', ...
    [1 0 0],'edgecolor', 'none')
plot3(orbit_Tot(:,1),orbit_Tot(:,2),orbit_Tot(:,3),'b-')
xlabel('distance (km)'),ylabel('distance (km)'),zlabel('distance (km)')

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
elapsed_time(i) = T*NumOrbs + time2(length(reentry));

end

DV = DV*V0 - V0;

figure(11) , hold on
leg = legend([num2str(DV(1))], [num2str(DV(2))],[num2str(DV(3))],...
        [num2str(DV(4))],[num2str(DV(5))],[num2str(DV(6))],...
        [num2str(DV(7))],[num2str(DV(8))],[num2str(DV(9))],...
        [num2str(DV(10))],'location','eastoutside');
title(leg,"deltaV (km/s)") 
xlabel('distance (km)'),ylabel('heating rate (kg/s^3)')

format short
t1 = table(DV(:),elapsed_time(:),vfinal(:));
t1.Properties.VariableNames = {'DeltaV', 'Duration', 'Final_Velocity'};

t2 = table(gamma_final(:),maxQ(:), alt_MaxQ(:));
t2.Properties.VariableNames = {'Flight_Path_Angle', 'Max_Heating_Rate',...
        'Altitude_MHR'};
    
figure(12)
uitable('Data',t1{:,:},'ColumnName',t1.Properties.VariableNames,...
        'Units','Normalized', 'Position',[0, 0, 1, 1],'FontSize',12);
    
figure(13)
uitable('Data',t2{:,:},'ColumnName',t2.Properties.VariableNames,...
        'Units','Normalized', 'Position',[0, 0, 1, 1],'FontSize',12);





