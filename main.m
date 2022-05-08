%% main.m
% @author: Agnello Hugo, Borbouse Maxime, Camus Sarah, Lucarelli Lorenzo, 
%          Miny Héloïse, Trifilò Tecla
% Date: May 2022
% Description: main code for the space experiment development project

%% Initialization
Orbit.rho = 0;
Orbit.R = 35000e3; % [km]
A = ones(4, 1);
minRequirements = [0 -40 -30 10 -30 ] + 273.15;
maxRequirements = [30 80 65 30 65] + 273.15;
Q = 0;
Mission.tau_eclipse = 1/3*24.6229*3600;
Mission.tau_sun = 2/3*24.6229*3600;
Mission.tau = 24.6229*3600;
powRequirements = [3.5 0.01 0.06 0.65 3.22; ...
                   0 0 0.06 6.5 3.22];


%% Run codes
[TCold, THot, alphaEpsilon] = ThermalDesign(Orbit, A, minRequirements, maxRequirements, Q, 1);
Array = PowerBudget(Mission, powRequirements);