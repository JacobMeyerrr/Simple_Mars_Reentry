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
%   USAGE:
%
% This function takes in an altitude and calculates the density in the
% Martian atmosphere to output for use in other functions. 
%
%   INPUTS: 
%                                                       UNITS:
%               h           Altitude                    [km]
%
%   OUTPUTS:
%
%               rho         Density found at altitude   [kg/m^3]
%
%   INPUT CHECKING:

if( ~isnumeric(altitude) )
    error('input altitude needs to be numeric!')
end
if( altitude < 0 )
        error('input altitude must be positive!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input checking to possible add? 
% if NaN is present?
% NEED TO DOUBLE CHECK THE UNITS

% Initialize h as the altitude value
h = altitude;

%% Calculate the Density on Mars and the altitude

% If altitude is less then 7km (7000m) (surface to mid-atmosphere)
if( 0 <= h ) && ( h <= 7)
    
    % Convert the altitude to meters
    h = altitude * 1000;

    % Calculate the pressure
    % P = pressure (kPa)
    P = 0.699*exp(-0.00009*h);
    
    % Calculate the temperature
    % T = temperature (deg C)
    T = -31-0.000998*h;
    
    % Calculate the lower density
    % NOTE: temp is in degree celsius but gets converted to kelvin here
    rho = P/(0.1921*(T+273.1));
    
%If altitude is greater than 7 km (7000 meters) and less than or equal to
%65000 meters (mid-atmosphere to upper atmosphere)
elseif( 7 < h ) && ( h <= 65 )
    
    % Convert the altitude to meters
    h = altitude * 1000;
    
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
elseif( 65 < h ) && ( h < 227920000 )
    
   % Declare coefficents on the polynomial density function
    a0 = 49.8118119899434;
    a1 = -5.9123700325916;
    a2 = -3.5638800977374;
    a3 = 0.380908561109888;
    
    % Calculate the upper density
    rho = 0.88325*exp(a0+a1*log(h)+a2*(log(h))^2+a3*(log(h))^3);
    
% If the orbit exceeds the semimajor axis of mars, according to
% NASA.GOV
elseif( h > 227920000 )
    warning("YOU HAVE ESCAPED THE SEMIMAJOR AXIS ORBIT OF MARS!!")
    
else(error( "An unknown error occured calculating or from the altitude input"))

end

% Output the density
density = rho;

