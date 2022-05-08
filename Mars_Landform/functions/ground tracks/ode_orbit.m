function dy = ode_orbit( ~, y, mu )
%ode_orbit ODE system for the orbit
%
% PROTOTYPE
% dy = ode_orbit( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[3x2] State of the s/c (position and velocity) [ L, L/T ]
% mu[1] Gravitational parameter [L^3/T^2]
%
% OUTPUT:
% dy[2x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% Set the derivatives of the state
t=sqrt(y(1)^2+y(2)^2+y(3)^2);
dy = [ y(4)
    y(5)
    y(6)
(-mu/t^3)*y(1)
(-mu/t^3)*y(2)
(-mu/t^3)*y(3)];
end