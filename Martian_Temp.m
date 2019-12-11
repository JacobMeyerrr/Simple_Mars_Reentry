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
%               can be a vector of heights, or a scalar value.
%
% OUTPUTS:
%       temp:   atmospheric temperature at given altitude  (degrees C)
%               will be a vector or a scalar depending on input

%% Function Main

h = height;

T = zeros(size(h));

h0_7  = find(h>0 & h< 7);
h7_65 = find(h>= 7 & h < 65);
h65p  = find(h > 65);


T(h65p)  = -167.7;
T(h7_65) = -23.4 - 0.00222 * h(h7_65);
T(h0_7)  = -31 - 0.000998 * h(h0_7);

temp = T;


