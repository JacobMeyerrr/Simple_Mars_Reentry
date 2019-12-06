function temp = Martian_Temp(height)
%% Martian_Temp
% calculates the temperature of the Martian atmosphere for a given
% altitude in meters.
%
% for height > 7 (km)
%       T = -23.4 - 0.00222*h*1000
%
% for height < 7 (km)
%       T = -31 - 0.000998*h*1000
%
% USAGE:
%       Martian_Temp(height)
%
% INPUTS:
%       height: altitude above msl of Mars  (km)
%
% OUTPUTS:
%       temp:   atmospheric temperature at given altitude  (degrees C)
%

%% Function Main

h = height * 1000;  % converts height from units of km to units of m.

if h >= 65000
    T = -167.7;
    
elseif h >= 7000
    T = -23.4 - 0.00222 * h;

else 
    T = -31 - 0.000998 * h;
end

temp = T;
