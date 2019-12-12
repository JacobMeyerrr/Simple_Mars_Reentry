
%% Mars Orbital Reentry Simulation
% Simulate the orbital path around Mars for the Nasa Orion capsual.
% Implement deorbit mavuever and simulate entry into the Martian atmosphere
% and terminate the simulation at 5km.
%
% Calculate the Body-Averaged heat rate during reentry
% Calculate the final velocity, flightpath angle, and elapsed time.
%

%% Defining Constants

%Planet constants
Radius_Mars = 3390;                     %[km]
Mass_Mars   = 6.39e23;                  %[kg]

%Spacecraft Constants
Mass_craft = 10400;                     %[kg]

%Gravitational Constant
G          = (6.67408e-11)/1000^3;      %[km^3/(kg*s^2)]

%Initial conditions
Altitude = Radius_Mars + 200;           %[km]
V0 = sqrt(G*(Mass_Mars + Mass_craft) ... 
              /Altitude)+0.0073062;     %[km/s]


%Calculating the orbital period for a circular orbit
T = sqrt(4*pi^2*(Altitude)^3 ...
                    /(G*Mass_Mars));
                
%Defining a time vector
NumOrbs = 3;                            %approximate number of orbits
time1 = linspace(0,T*NumOrbs,5e5);      %[s]

%% Calculating deltaV required to induce reentry
% Define a desired periapsis that will induce a reentry.
% run simulation with varying deltaV's and a zero drag force model
% measuring the periapsis.
% exit simulation when the desired periapsis is achieved
%

dv = linspace(1,0.5,250);                 %Vector defining percent V0
entry = 50;                               %Desired periapsis [km]
enters = 0;

for i = 1:length(dv)

    V = V0*dv(i);

    vHat     = [0;1/sqrt(2);1/sqrt(2)];
    PosState = [Altitude;0;0];            %[x;y;z]   [km;km;km]
    VelState = V*vHat;                    %[Vx;Vy;Vz][km/s;km/s;km/s]

    state0 = [PosState; VelState];        %initial state vector

    %Numerical Integration - calculating the orbit after deorbit burn.
    orbit = ODENumIntRK4(@CraftOrbit,time1,state0,Mass_craft,false);

    R = vecnorm(orbit(:,1:3),2,2)-Radius_Mars;

    periapsis = min(R);
    
    %If the periapsis of the previous iteration meets the desired value -
    %print the altitude vs. time to verify and break the loop.
    if(periapsis <= entry)
        x = ones(size(R))*entry;
        figure(1), hold on
        plot(time1,x,'r-'),plot(time1,R,'b-')
        title("Orbital Altitude vs. Time")
        xlabel("time (s)"),ylabel("Altitude (km)")
        legend("Target Periapsis","Orbital Altitude")
        initialDV = dv(i);
        break;
    end
    
    %If the periapsis of the previous iteration was within 20% of the
    %desired value redefine dv with smaller divisions to increase 
    %accuracy.
    if(periapsis <= entry*1.2 && enters == 0)
        dv = linspace(dv(i),dv(i)/2,1e4);
        i=1;
        enters = 1;
    end
end


%% Simulate reentry into Martian atmosphere

%Initial conditions
Altitude = Radius_Mars + 200;           %[km]
V0 = sqrt(G*(Mass_Mars+Mass_craft) ... 
              /Altitude)+0.0073062;     %[km/s]

%Define initial state vector and time vector
vHat     = [0;1/sqrt(2);1/sqrt(2)];     %initial velocity unit vector
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

% Calculate initial parking orbit using numerical integration
orbit = ODENumIntRK4(@CraftOrbit,time1,state0,Mass_craft,true);

%Initialize deorbit burn
v = orbit(end,4:6);             % current velocity state
u = norm(v);                    % magnitude of velocity
vHat = v/u;                     % velocity unit vector

deltaV = DV(i);                 % DeltaV to induce reentry

Vnew = u*deltaV*vHat;           % New Velocity vector after DeltaV impulse

orbit(end,4:6) = Vnew;          % Orbit(end,:) is our new state vector

time2 = linspace(0,1e6,1e6);    % Define new time vector

%Simulate spacecraft reentry into Mars atmosphere
reentry = ODENumIntRK4(@CraftOrbit,time2,orbit(end,:),Mass_craft,true);

orbit_Tot = [orbit;reentry];    % Matrix of states for entire duration

%Plotting the orbital path around Mars
figure(i+1), hold on, grid on
[x,y,z] = sphere(50);
surf(x*Radius_Mars,y*Radius_Mars,z*Radius_Mars,'FaceColor', ...
    [1 0 0],'edgecolor', 'none')
plot3(orbit_Tot(:,1),orbit_Tot(:,2),orbit_Tot(:,3),'b-')
title("Orbital Trajectory around Mars")
xlabel('distance (km)'),ylabel('distance (km)'),zlabel('distance (km)')

%Calculate Body-Averaged HeatingRate
qAvg = HeatingRate(orbit_Tot);

R = vecnorm(orbit_Tot(:,1:3),2,2)-Radius_Mars;

%Plotting every iterations Heating rate vs. Altitude
figure(12), grid on, hold on
plot(R,qAvg)

%Calculate Max Heating Rate, and Alt where it occurs
maxQ(i) = max(qAvg);
alt_MaxQ(i) = R(find(qAvg == maxQ(i)));

%Record Final Velocity, Flight Path Angle (Gamma), and time elapsed
vfinal(i) = norm(orbit_Tot(end,4:6));
gamma_final(i) = flightpathangle(orbit_Tot(end,1:3),orbit_Tot(end,4:6));
elapsed_time(i) = T*NumOrbs + time2(length(reentry))/3600;

end

% Calculating the deltaV used to induce reentry
DV = DV*V0 - V0;

%Formatting the figure for Heating rate vs. Altitude. 
figure(12) , hold on
leg = legend([num2str(DV(1))], [num2str(DV(2))],[num2str(DV(3))],...
        [num2str(DV(4))],[num2str(DV(5))],[num2str(DV(6))],...
        [num2str(DV(7))],[num2str(DV(8))],[num2str(DV(9))],...
        [num2str(DV(10))],'location','eastoutside');
title(leg,"deltaV (km/s)") 
title("Heating Rate vs. Altitude")
xlabel('Altitude (km)'),ylabel('heating rate (kg/s^3)')


%Building tables of final data to display
format short
t1 = table(DV(:),elapsed_time(:),vfinal(:));
t1.Properties.VariableNames = {'DeltaV', 'Duration', 'Final_Velocity'};

t2 = table(gamma_final(:),maxQ(:), alt_MaxQ(:));
t2.Properties.VariableNames = {'Flight_Path_Angle', 'Max_Heating_Rate',...
        'Altitude_MHR'};
    
figure(13)
uitable('Data',t1{:,:},'ColumnName',t1.Properties.VariableNames,...
        'Units','Normalized', 'Position',[0, 0, 1, 1],'FontSize',13);
    
figure(14)
uitable('Data',t2{:,:},'ColumnName',t2.Properties.VariableNames,...
        'Units','Normalized', 'Position',[0, 0, 1, 1],'FontSize',13);

%Plotting max heating rate vs. deltaV used.
figure(15)
plot(DV,maxQ,'r.','markersize',10)
title("Max Heating Rate vs. deltaV")
xlabel("deltaV (km/s)"),ylabel("Max Heating Rate (kg/s^3)")
grid on, grid minor



