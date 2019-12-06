function s = Mars_SpeedofSound(temp)
%% Mars_SpeedofSound
%
% Calculates the speed of sound in martian atmosphere given the
% temperature.
%
% Assumes atmosphere composition of 100% CO2
%
% USAGE:
%       Mars_SpeedofSound(temp)
%       
% CONSTANTS:
%       gamma = 1.28      -  adiabatic constant for CO2  (unitless)  
%           R = 8.314462  -  Universal gas constant      (kg*m/K*mol*s^2) 
%           M = 44.01e-3  -  Molecular Mass of CO2       (kg/mol)
%
% INPUTS:
%       temp:   current temperature of martian atmosphere
%
% OUTPUTS:
%       s:      speed of sound in martian atmosphere
%
% 
%% Function Main

%Define constants
gamma = 1.28;       %unitless
R = 8.314462;       %kg*m^2/k*mol*s^2 
M =  44.01e-3;      %(kg/mol)

%Convert input temperature to from celsius to Kelvin
T = temp + 273.15;  %K

%Calculate Speed of Sound
s = sqrt(gamma*R*T/M);
