% Accuracy computation
e = .7; % [km]
FoV = 9; % [deg]
h = 500; % [km]
s = 99.3; % [km]

e = h*tan(theta + FoV/2)-tan(FoV/2);
theta = atan((s/2 + e)/h)-atan(s/(2*h));
theta = rad2deg(theta);

% disturbance torques
J = 590;
c = 3e8;
As = 0.06;
q = 0.3*0.2+0.3*0.1+0.2*0.1;
DeltaX = 0.01;
muM = 42828.37;
R = 3389.5 + 500;
I1 = 0.3^2+0.2^2; % with a 12kg 6u cubesat
I2 = 0.3^2+0.1^2;
I3 = 0.2^2+0.1^2;
DeltaI = I1-I3;
Pgg = 2*pi*sqrt(R^3/muM);
Psrp = 0.8*Pgg;

T_gg = 3*muM/2/R^3*DeltaI;
T_SRP = J/c*As*(1+q)*DeltaX;

H_gg = 0.707*T_gg*Pgg/4;
H_SRP = 0.707*T_SRP*Psrp/4;