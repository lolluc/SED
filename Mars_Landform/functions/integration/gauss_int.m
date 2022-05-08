function [KEP,t]= gauss_int(KEP_0,t,A_M,Cd,ode_num,mu)
% 
% Propagate the orbit using gauss planetary equations. Only j2 and drag
% disturbances are modeled.
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
%   KEP[Nx6]      history of propagated keplerian elements [a,e,i,OM,om,th]
%                  
%   t  [1xN]      time vector associated to KEP [s]
%
% REQUIRED FUNCTIONS: astroConstants, rho_atm
%  
%

if nargin==4
    ode_num='113';
    mu=astroConstants(13);
elseif nargin==5
    mu=astroConstants(13);
end


J2=0.00108263;
R_e=6378.137;
omega_e=0.000072921;

options= odeset('reltol',1e-13,'AbsTol',1e-14);

switch ode_num %chose ode solver
    
    case '45'
        [t,KEP]= ode45(@(~,K) ode_gauss(t,K),t,KEP_0,options);
        
    case '113'
        [t,KEP]= ode113(@(~,K)ode_gauss(t,K),t,KEP_0,options);
end

% ODE function
function [dkepdt]=ode_gauss(~,K)
a=K(1); e=K(2); i=K(3); OM=K(4); om=K(5); th=K(6);  u=om+th;

% compute useful parameters to reduce computations
p=a*(1-e^2);
h= sqrt(mu*p);
r= p/(1+e*cos(th));

%j2 accelerations in rsw frame
k_j2=-3/2*J2*mu*R_e^2/r^4;
a_j2_r= k_j2*(1-3*(sin(i)^2)*(sin(u)^2));
a_j2_s= k_j2*(sin(i)^2*sin(2*(u)));
a_j2_w= k_j2*(sin(2*i)*sin(u));

% kep2car algorithm to get velocity and position in cartesian coordinates
r_pf = r*[cos(th); sin(th); 0];
v_pf = sqrt(mu/p)*[-sin(th); e + cos(th); 0];
R3_OM = [cos(OM) sin(OM) 0 ; -sin(OM) cos(OM) 0 ; 0 0 1];
R1_i = [1 0 0 ; 0 cos(i) sin(i) ; 0 -sin(i) cos(i)];
R3_om = [cos(om) sin(om) 0 ; -sin(om) cos(om) 0 ; 0 0 1];
T_mat = (R3_om*R1_i*R3_OM)'; 
r_vect = T_mat*r_pf;
v_vect = T_mat*v_pf;

%relative velocity wrt the atmosphere and drag vector
v_rel=[v_vect(1)+omega_e*r_vect(2);v_vect(2)-omega_e*r_vect(1);v_vect(3)];
a_drag_xyz=-0.5*1e3*rho_atm(r-R_e)*A_M*Cd*norm(v_rel)*v_rel;

%construction of rsw frame
r_dir= r_vect/norm(r_vect);
w_dir= cross(r_dir,v_vect);w_dir=w_dir/norm(w_dir);
s_dir=cross(w_dir,r_dir);

%rotation of drag vector in rsw frame
a_drag_rsw=[r_dir';s_dir';w_dir']*a_drag_xyz;

%total accelerations in rsw frame
a_r=a_drag_rsw(1)+a_j2_r;
a_s= a_drag_rsw(2)+a_j2_s;
a_w= a_drag_rsw(3)+ a_j2_w;

%gauss planetary equations in rsw frame
a_dot= 2*a^2/h*(e*sin(th)*a_r+p/r*a_s);
e_dot= 1/h*(p*sin(th)*a_r+((p+r)*cos(th)+r*e)*a_s);
i_dot= r/h*cos(u)*a_w;
OM_dot=r/h*sin(u)/sin(i)*a_w;
om_dot=1/(h*e)*(-p*cos(th)*a_r+(p+r)*sin(th)*a_s)-r/h*sin(u)/tan(i)*a_w;
th_dot= h/r^2+1/(e*h)*(p*cos(th)*a_r-(p+r)*sin(th)*a_s);

dkepdt= [a_dot;e_dot;i_dot;OM_dot;om_dot;th_dot];

end
end
        