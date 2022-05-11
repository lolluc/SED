function [r,v,t]= cart_int(KEP_0,t,A_M,Cd,ode_num,mu)
% 
% Propagate the orbit using kepler planetary equations in cartesian 
% coordinates ì. Only j2 and drag disturbances are modeled.
%
%
% INPUTS:
%   KEP_0 [1x6]   keplerian elements of the initial condition, where
%                 KEP_0=[a,e,i,OM,om,th] [km,-,rad,rad,rad,rad]
% 
%   t [1xN]       time integration vector. Can be [1x2] and include only
%                 initial and final time or [1xN] to chose wanted time 
%                 instants. [s] 
%
%   A_M[1x1]      Area to mass ratio of the satellite. [m^2/kg]
%
%   Cd [1x1]      Drag coefficient of the satellite. [-]
%
%   ode_num       string with code of the wanted solver: (default '113')
%                 -'45'  to use ode45
%                 -'113' to use ode113  
%  
%   mu [1x1]      central body planetary constant (default earth)[km^3/s^2]
% 
% 
% OUTPUTS:
%   r[Nx3]        history of propagated radius vectors [km,km,km]
%
%   v[Nx3]        history of propagated velocity vectors [km/s,km/s,km/s]
%                  
%   t  [1xN]      time vector associated to r and v [s]
%
% REQUIRED FUNCTIONS: astroConstants, rho_atm, kep2car
%  
%
if nargin==4
    ode_num='113';
    mu=astroConstants(13);
elseif nargin==5
    mu=astroConstants(13);
end


%gather intial condition state vector with kep2car
[r0,v0]= kep2car(KEP_0(1),KEP_0(2),KEP_0(3),KEP_0(4),KEP_0(5),KEP_0(6),mu);


J2=0.00108263;
R_e=6378.137;
omega_e=0.000072921;

options= odeset('reltol',1e-13,'AbsTol',1e-14);

switch ode_num %select ode solver
    
    case '45'
        [t,u]= ode45( @(~,r) ode_cart(t,r), t,[r0;v0],options);
        r= u(:,1:3);
        v= u(:,4:6);
    case '113'
        [t,u]= ode113( @(~,r) ode_cart(t,r), t,[r0;v0],options);
        r= u(:,1:3);
        v= u(:,4:6);
end

%ODE
function [drdt]=ode_cart(~,r)
%compute radius to save time
R= sqrt(r(1)^2+r(2)^2+r(3)^2);   

%compute j2 acceleration
j2_con=3/2*J2*mu*R_e^2/(R^4);
j2x=j2_con*(r(1)/R*(5*r(3)^2/R^2-1));
j2y=j2_con*(r(2)/R*(5*r(3)^2/R^2-1));
j2z=j2_con*(r(3)/R*(5*r(3)^2/R^2-3));

%compute relative velocity wrt the atmosphere and drag acceleration
v_rel=[r(4)+omega_e*r(2);r(5)-omega_e*r(1);r(6)];
a_drag= -0.5*1e3*rho_atm(R-R_e)*A_M*Cd*norm(v_rel)*v_rel;

%output state derivatives
drdt=[r(4); r(5) ; r(6); ... 
    -mu*r(1)/(R^(3))+j2x+a_drag(1);... 
    -mu*r(2)/(R^(3))+j2y+a_drag(2);... 
    -mu*r(3)/(R^(3))+j2z+a_drag(3)];
end
end






