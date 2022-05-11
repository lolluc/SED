% Accuracy computation
e = 2; % [km]
FoV = ;
h = ;

e = h*tan(theta + FoV/2)-tan(FoV/2);
theta = atan((s/2 + e)/h)-atan(s/(2*h));

T_gg = 3*muM/2/R^3*DeltaI;
T_SRP = J/c*As*(1+q)*DeltaX;

H_gg = 0.707*T_gg*P/4;
H_SRP = 0.707*T_SRP*P/4;