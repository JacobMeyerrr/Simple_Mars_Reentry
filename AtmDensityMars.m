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
        fprintf("altitude:%d",altitude)
        error('input altitude must be positive');
end
% input checking to possible add? 
% if NaN is present, if there is a negative altitude?

% Initialize h as the altitude value
h = altitude;
h0_7  = find(h>0 & h<=7);
h7_65 = find(h>7 & h<=65);
h65p  = find(h>65);

% Initalize density vector
rho = zeros(size(h));
P = zeros(size(h));
%% Calculate the Density on Mars and the altitude

if size(h0_7) >= [1,1] % Checks if there's any altitudes between the surface and 7 km
% If altitude is less then 7km (7000m) (surface to mid-atmosphere)
    
    % Convert the altitude to meters
    h(h0_7) = h(h0_7) * 1000;

    % Calculate the pressure
    % P = pressure (kPa)
    P(h0_7) = 0.699*exp(-0.00009*h(h0_7));
    
    % Calculate the temperature
    % T = temperature (deg C)
    T(h0_7) = Martian_Temp(h(h0_7));
     
    T = T'; % Convert row to column vector
    % Calculate the lower density
    % NOTE: temp is in degree celsius but gets converted to kelvin here
    rho(h0_7) = P(h0_7)./(0.1921*(T(h0_7)+273.1));
end

if size(h7_65) >= [1,1] % Checking if there's any altitudes above 7 km, below 65 km
%If altitude is greater than 7 km (7000 meters) and less than or equal to
%65000 meters (mid-atmosphere to upper atmosphere)

    % Convert the altitude to meters
    h(h7_65) = h(h7_65) * 1000;
    
    % Calculate the pressure
    % P = pressure (kPa)
    P(h7_65) = 0.699*exp(-0.00009*h(h7_65));
    
    % Calculate the temperature
    % T = temperature (deg C)
    T(h7_65) = Martian_Temp(h(h7_65));
    
    % Calculate the lower density
    % NOTE: temp is in degree celsius but gets converted to kelvin here
    rho(h7_65) = P(h7_65)./(0.1921*(T(h7_65)+273.1));
end 

if size(h65p) >= [1,1] % If there's any altitudes above 65 km
% If altitude is greater than 65000 meters (upper atmosphere to space)
    
   % Declare coefficents on the polynomial density function
    a0 = 49.8118119899434;
    a1 = -5.9123700325916;
    a2 = -3.5638800977374;
    a3 = 0.380908561109888;
    
    % Calculate the upper density
    rho(h65p) = 0.88325*exp(a0+a1*log(h(h65p))+a2*(log(h(h65p))).^2 ...
                        +a3*(log(h(h65p))).^3);
end    
% else
%     disp(altitude)
%     error( "An unknown error occured calculating or from the altitude input")

% Output the density

density = rho;

