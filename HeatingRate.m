function qAvg = HeatingRate(States)

%HeatingRate
%
%Usage:
%   This function is used to calculate the heating rate of the capsule
%   during the reentry period.
%
%   qAvg = HeatingRate(States) 
%
%Input:
%   States      A vector containing the position and velocity values of the
%               capluse.
%
%Output:
%   qAvg        The average heating rate during reentry     [kg/s^3]
%

%noramize radius vector
R = vecnorm(States(:,1:3),2,2);

%normalize velocity vector
V = vecnorm(States(:,4:6),2,2);

%calculate atmospheric density based on altitude
rho = AtmDensityMars(R-3390);

%calculate temperature based on altitude
Temp = Martian_Temp(R);

%calculate kinematic viscosity based on temperature and atmospheric
%density
nu = KinVisc(Temp,rho);

%calculate  speed of sound based on current temperature
S = Mars_SpeedofSound(Temp);

%calculate  mach number based on velocity and speed of sound
M = Mach_Number(V,S);

%calculate the Reynolds number based on velocity, radius of capsule, and
%kinematic viscosity
Re = Reynolds_Number(V,5,nu);

%calculate the body average skin friction coefficient based on the Reynolds
%number and the Mach number
Cf = BA_SkinFric(Re,M);

%calculate the average heating rate based on the atmospheric density,
%velocity, and the body avreage skin friction coefficient.
qAvg = BA_HeatRate(rho,V,Cf);
