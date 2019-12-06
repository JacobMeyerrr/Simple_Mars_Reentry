function nu = KinVisc(temp,rho)

%KinVisc
%Usage:
%   KinVisc computes the kinematic viscosity of the Martian atmosphere at a
%   specific temperature and atmospheric density. 
%
%Input:
%   temp: temperature of the atmosphere    [oC]
%   rho: atmospheric density               [kg/m^3]
%
%Output:
%   v: kinematic viscosity                 [m^2/s]
%

%constants
B=1.458*10^-6;
S=110.4;
T=temp+273.15;

%dynamic viscosity
u=(B*T^1.5)/(T+S);

nu=u/rho;
