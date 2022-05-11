%% main.m
% @author: Agnello Hugo, Borbouse Maxime, Camus Sarah, Lucarelli Lorenzo, 
%          Miny Héloïse, Trifilò Tecla
% Date: May 2022
% Description: main code for the space experiment development project

%% Initialization
Orbit.rho = 0;
Orbit.R = 35000e3; % [km]
A = [0.3*0.2 0.3*0.2 0.3*0.2 0.3*0.2];
minRequirements = [10 -40 -30 10 -30] + 273.15; % 1. Optic, 2. 
maxRequirements = [30 80 65 30 65] + 273.15;
Q = 0;
Mission.tau_eclipse = 3500;
Mission.tau_sun = 8800 - 3500;
Mission.tau = 8800;
Mission.m = 12;
powRequirements = [4.5 0.01 0.06 0.65 3.22; 0 0 0.06 6.5 3.22]; % 1. Optic, 2.

%% Run codes
[TCold, THot] = ThermalDesign(Orbit, A, 0.26, 0.89, minRequirements, maxRequirements, Q, 0);
Array = PowerBudget(Mission, powRequirements);