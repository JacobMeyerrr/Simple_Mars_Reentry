function density = AtmDensityMars(altitude)
%========================================================================%
%   m-script name:        AtmDensityMars.m
%   Programmer:           Dan Bombeck
%   Date Created:         12-6-2019   
%   Date Last Modified:   
%   
%   Description:          This m-script finds density on Mars using
%                         altitude
%========================================================================%
%
% This function takes in an altitude and calculates the density in the
% Martian atmosphere to output for use in other functions. It evaluates 
%
%   INPUTS: 
%                                                       UNITS:
%               h           Altitude                    [Km]
%
%   OUTPUTS:
%
%               rho         Density found at altitude   [Kg/m^3]
%
%   INPUT CHECKING:

if( ~isnumeric(altitude) )
    error('input altitude needs to be numeric!')
end
% input checking to possible add? 
% if NaN is present, if there is a negative altitude?

%% Convert the altitude to meters
h = altitude * 1000;

%% Calculate the Density on Mars and the altitude

% If altitude is less then 7000 meters (surface to mid-atmosphere)
if( 0 <= h ) && ( h <= 7000)
    
    % Calculate the pressure
    % P = pressure (kPa)
    P = 0.699*exp(-0.00009*h);
    
    % Calculate the temperature
    % T = temperature (deg C)
    T = -31-0.000998*h;
    
    % Calculate the lower density
    % NOTE: temp is in degree celsius but gets converted to kelvin here
    rho = P/(0.1921*(T+273.1));
    
%If altitude is greater than 7000 meters and less than or equal to
%65000 meters (mid-atmosphere to upper atmosphere)
elseif( 7000 < h ) && ( h <= 65000 )
    
    % Calculate the pressure
    % P = pressure (kPa)
    P = 0.699*exp(-0.00009*h);
    
    % Calculate the temperature
    % T = temperature (deg C)
    T = -23.4-.00222*h;
    
    % Calculate the lower density
    % NOTE: temp is in degree celsius but gets converted to kelvin here
    rho = P/(0.1921*(T+273.1));
    
% If altitude is greater than 65000 meters (upper atmosphere to space)
elseif( 65000 < h )
    
   % Declare coefficents on the polynomial density function
    a0 = 49.8118119899434;
    a1 = -5.9123700325916;
    a2 = -3.5638800977374;
    a3 = 0.380908561109888;
    
    % Calculate the upper density
    rho = 0.88325*exp(a0+a1*log(h)+a2*(log(h))^2+a3*(log(h))^3);
    
else( "An unknown error occured")
end

% Output the density
density = rho;