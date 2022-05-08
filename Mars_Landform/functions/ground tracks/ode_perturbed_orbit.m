function dy = ode_perturbed_orbit( ~, y, mu, Re, j2 )
%ode_perturbed_orbit ODE system for the perturbed 2BP
%
% PROTOTYPE
% dy = ode_perturbed_orbit( t, y, mu, Re, j2 )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[3x2] State of the s/c (position and velocity) [ L, L/T ]
% mu[1] Gravitational parameter [L^3/T^2]
% Re[1] Radius of Earth [L]
% j2[1] Second zonal armonic [-]
%
% OUTPUT:
% dy[2x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% Set the derivatives of the state
t=sqrt(y(1)^2+y(2)^2+y(3)^2);
a1=3/2*j2*mu*Re^2/t^4*y(1)/t*(5*y(3)^2/t^2-1);
a2=3/2*j2*mu*Re^2/t^4*y(2)/t*(5*y(3)^2/t^2-1);
a3=3/2*j2*mu*Re^2/t^4*y(3)/t*(5*y(3)^2/t^2-3);
dy = [ y(4)
    y(5)
    y(6)
(-mu/t^3)*y(1)+a1
(-mu/t^3)*y(2)+a2
(-mu/t^3)*y(3)+a3];
end