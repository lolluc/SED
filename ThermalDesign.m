%% ThermalDesign.m
% @author: Maxime Borbouse
% Date: May 2022
% Description: computes the temperature of the satellite in both cold and
%   hot cases by performing a power balance (1 node).
% Inputs: * Orbit: structure containing orbit parameters
%         * A [m^2] : array containing surface areas (1. projected area 
%              receiving solar radiation, 2. projected area receiving 
%              albedo radiation, 3. projected area receiving planetary 
%              radiation, 4. total area)
%         * minRequirements: array containing minimum temperature
%                            requirements for each component
%         * maxRequirements: array containing maximum temperature
%                            requirements for each component
%         * showplot: boolean (1 to show the plots, 0 to hide them)
% Outputs: * TCold [K]: Satellite temperature in cold case
%          * THot [K]: Satellite temperature in hot case

function [TCold, THot, alphaEpsilon] = ThermalDesign(Orbit, A, minRequirements, maxRequirements, Q, showplot)
%% THERMAL ENVIRONMENT
% Solar radiation
Js = parameters.P/(4*pi*(parameters.D*1000)^2); % [W/m^2] Solar radiation intensity

% Albedo radiation
d = Orbit.R/parameters.R;
x = sqrt(d^2 - 1);
y = -x*cot(Orbit.rho);
if abs(Orbit.rho) < pi/2 - asin(1/d)
    F = cos(Orbit.rho)/d^2; % [-] Visibility factor
else
    F = 1/(pi*d^2)*(cos(Orbit.rho)*acos(y) - c*sin(Orbit.rho*sqrt(1-y^2)) + 1/pi*atan(sin(Orbit.rho)*sqrt(1 - y^2)/x));
end
Ja = Js*parameters.albedo*F; % Intensity of albedo radiation

% Planetary radiation
J_M = parameters.eps*parameters.sigma*parameters.Tb; % [W/m^2] Mars thermal radiation
Qir = J_M;
Rrad = parameters.R; % [m] Radius of Mars effective radiating surface
Jp = Qir*(Rrad/Orbit.R)^2;

%% THERMAL BALANCE
Asolar = A(1); % [m^2] Projected area receiving solar radiation
Aalbedo = A(2); % [m^2] Projected area receiving albedo radiation
Aplanetary = A(3); % [m^2] Projected area receiving planetary radiation
Asurface = A(4); % [m^2] Total area

for alphaEpsilon = 0:0.001:1
    if Q ~= 0
        for epsilon = 0.001:0.001:1
            TCold = (Aplanetary*Jp/(Asurface*parameters.sigma) + Q/(Asurface*parameters.sigma*epsilon) + ...
                (Asolar*Js + Aalbedo*Ja)/(Asurface*parameters.sigma)*(alphaEpsilon)).^(1/4); % [K] Spacecraft temperature
            THot = (Aplanetary*Jp/(Asurface*parameters.sigma) + Q/(Asurface*parameters.sigma*epsilon) + ...
                Aalbedo*Ja/(Asurface*parameters.sigma)*(alphaEpsilon)).^(1/4); % [K] Spacecraft temperature
        end
    else
        TCold = (Aplanetary*Jp/(Asurface*parameters.sigma) + ...
        (Asolar*Js + Aalbedo*Ja)/(Asurface*parameters.sigma)*(alphaEpsilon)).^(1/4); % [K] Spacecraft temperature
    THot = (Aplanetary*Jp/(Asurface*parameters.sigma) + ...
        Aalbedo*Ja/(Asurface*parameters.sigma)*(alphaEpsilon)).^(1/4); % [K] Spacecraft temperature
    end
    
    if TCold > mean(minRequirements) && THot < mean(maxRequirements)
        fprintf(['<strong> Absorptivity/Emissivity:</strong>           alpha/epsilon = ' num2str(alphaEpsilon) ' [-]\n'])
        fprintf(['<strong> Cold case temperature:</strong>             TCold = ' num2str(TCold) ' K\n'])
        fprintf(['<strong> Hot case temperature:</strong>              THot = ' num2str(THot) ' K\n'])
        break;
    end
    if TCold > mean(minRequirements) && THot < mean(maxRequirements)
        break;
    end
end

% Temperature with respect to absorptance and emissivity
if showplot
    alphas = 0:0.001:1;
    epsilons = 1;
    T_ae = (Aplanetary*Jp/(Asurface*parameters.sigma) + Q./(Asurface*parameters.sigma*epsilons) + ...
        (Asolar*Js + Aalbedo*Ja)/(Asurface*parameters.sigma)*(alphas./epsilons)).^(1/4); % [K] Spacecraft temperature

    figure(1)
    plot(alphas, T_ae, 'LineWidth', 1);
    xlabel('$\alpha = \epsilon$ [-]', 'Interpreter', 'Latex')
    ylabel('$T$ [K]', 'Interpreter', 'Latex')
end
end