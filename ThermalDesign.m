function T = ThermalDesign(Orbit, A, minRequirements, maxRequirements, Q, showplot)
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
for alpha = 0:0.001:1
    for epsilon = 0:0.001:1
        TCold = (Aplanetary*Jp/(Asurface*sigma) + Q/(Asurface*sigma*epsilon) + ...
            (Asolar*Js + Aalbedo*Ja)/(Asurface*sigma)*(alpha/epsilon)).^(1/4); % [K] Spacecraft temperature
        THot = (Aplanetary*Jp/(Asurface*sigma) + Q/(Asurface*sigma*epsilon) + ...
            Aalbedo*Ja/(Asurface*sigma)*(alpha/epsilon)).^(1/4); % [K] Spacecraft temperature
        
        if TCold > max(minRequirements) && TCold < min(maxRequirements) && THot > max(minRequirements) && THot < min(maxRequirements)
            fprintf(['<strong> Absorptivity:     alpha = ' num2str(alpha) ' [-]\n'])
            fprintf(['<strong> Emissivity:       epsilon = ' num2str(emissivity) ' [-]\n'])
            fprintf(['<strong> Temperature:      T = ' num2str(T) ' K\n'])
            break;
        end
    end
end

% Temperature with respect to absorptance and emissivity
if showplot
    alphas = 0:0.001:1;
    epsilons = 1;
    T_ae = (Aplanetary*Jp/(Asurface*sigma) + Q./(Asurface*sigma*epsilons) + ...
        (Asolar*Js + Aalbedo*Ja)/(Asurface*sigma)*(alphas./epsilons)).^(1/4); % [K] Spacecraft temperature

    figure(1)
    plot(alphas, T_ae, 'LineWidth', 1);
    xlabel('$\alpha = \epsilon$ [-]', 'Interpreter', 'Latex')
    ylabel('$T$ [K]', 'Interpreter', 'Latex')
end
end