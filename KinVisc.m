function nu = KinVisc(T,rho)

%KinVisc
%Usage:
%   KinVisc computes the kinematic viscosity of the Martian atmosphere at a
%   specific temperature and atmospheric density. 
%
%Input:
%   T: temperature of the atmosphere    [K]
%   rho: atmospheric density            [kg/m^3]
%
%Output:
%   v: kinematic viscosity              [m^2/s]
%

%constants
B=1.458*10^-6;      % [kg/(s*m*K^.5)]
S=110.4;            % [K]

%Convert degrees Celsius to Kelvin
T = T + 273.15;     % [K]

%dynamic viscosity
u=(B*T.^1.5)./(T+S);

nu=u./rho;
