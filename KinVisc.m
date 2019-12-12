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

%constants in Sutherland equation
B=1.458*10^-6;      % [kg/(s*m*K^.5)]
S=110.4;            % [K]

%Convert degrees Celsius to Kelvin
T = T + 273.15;     % [K]

%dynamic viscosity calculated via the Sutherland equation
u=(B*T.^1.5)./(T+S);

%calculate kinematic viscosity
nu=u./rho;
