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

%% Run codes
[TCold, THot] = ThermalDesign(Orbit, A, minRequirements, maxRequirements, Q, 0);