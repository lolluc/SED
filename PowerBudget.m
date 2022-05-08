function Array = PowerBudget(Mission, Requirements)
% tbu
etaBCR = 0.9;
etaBDR = 0.9;
etaAR = 0.9;
etaCharge = 0.9;
etad = 0.9;
DOD = 0.1; % Depth of discharge
D = 0.1; % Degradation factor (see 10.3.1)
etaCell = 0.9; % Solar cell efficiency

deltaTheta = 1; % [°] Array pointing error with respect to the sun
etaPacking = 0.9; % [-] Cell packing efficiency

S = parameters.P/(4*pi*(parameters.D*1000)^2); % [W/m^2] Solar radiation intensity

eta = etaBCR*etaBDR*etaAR;

PCharge = 1/eta*sum(Requirements(:,1))*Mission.tau_eclipse/Mission.tau_sun;
Array.P = sum(Requirements(:,2)) + PCharge;

Array.E_B = sum(Requirements(:,1))*(Mission.tau - Mission.tau_sun)/(etaCharge*DOD);
Array.epsilon = Array.P*Mission.tau;
Array.A = Array.P/(S*cosd(deltaTheta)*etaCell*etaPacking*(1-D));
end