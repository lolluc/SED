function [alpha, delta, lon, lat] = perturbed_groundTrack(kep_el, tspan, mu, omega_m, t0, th_G0, R_mars, j2)
% perturbed_groundTrack computes the ground track of an orbit with the j2
% effect
%
% PROTOTYPE
% [alpha, delta, lon, lat] = perturbed_groundTrack(kep_el, lon_0, tspan, omegaE, t0, th_G0)
%
% INPUT:
% kep_el   [1x6] State of the orbit at initial time (keplerian elements)            [km, -, rad]
% lon_0    [1]   Longitude of Greenwich meridian at initial time                    [rad]
% tspan    [:x1] Vector of times at which the ground track will be computed         [s]
% mu       [1]   Gravitational parameter                                            [km^3/s^2]
% omegaE   [1]   Planet's rotation velocity                                         [deg/h]
% t0       [1]   Reference initial time                                             [s]
% th_G0    [1]   Initial position of the rotating planet with respect to x axis     [rad]
% Re       [1]   Radius of the planet                                               [km]
% j2       [1]   Perturbation parameter                                             [-]
%
% OUTPUT:
% alpha    [:,1] Right ascension in Earth Centered Equatorial Inertial frame        [rad]
% delta    [:,1] Declination in Earth Centered Equatorial Inertial frame            [rad]
% lon      [:,1] Longitude with respect to rotating Earth (0 at Greenwich meridian) [rad]
% lat      [:,1] Latitude with respect to rotating Earth                            [rad]
%
% Gives the ground track of an orbit perturbed with j2 effect

[r0,v0]=kep2car(kep_el(1),kep_el(2),kep_el(3),kep_el(4),kep_el(5),kep_el(6),mu);

% Initial conditions
y0 = [r0'; v0'];

% Set options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);

% Perform the integration
[t, Y1] = ode113(@(t, y) ode_perturbed_orbit(t,y,mu,R_mars,j2), tspan, y0, options);

% Find state vectors
r_orbit = [Y1(:,1),Y1(:,2),Y1(:,3)];



%pre allocation
delta=zeros(length(t),1);
alpha=delta;
th_gw=delta;
lon=delta;


x=r_orbit(:,1);
y=r_orbit(:,2);
z=r_orbit(:,3);
R = sqrt(x.^2+y.^2+z.^2);
%delta=asin(z./R);
%lat=wrapToPi(lat);

% alpha=zeros(length(tspan),1);
% for j=1:length(tspan)
%     if y(j)/R(j)>0
%         alpha(j)=acos(x(j)/R(j)/cos(delta(j)));
%     else
%         alpha(j)=2*pi-acos(x(j)/R(j)/cos(delta(j)));
%     end
% end
omega_m=deg2rad(omega_m/3600); % mars's rotational speed


for ii=1:length(t)
    
    delta(ii)= asin(r_orbit(ii,3)/norm(r_orbit(ii,:))); %declination
    
    alpha(ii) =atan2(r_orbit(ii,2),r_orbit(ii,1)); % rigth ascension
    
    th_gw(ii)=th_G0+omega_m*(t(ii)); % teta of greenwich (time)
    th_gw(ii)=wrapTo2Pi(th_gw(ii));
    
    lon(ii)= alpha(ii)-th_gw(ii); %longitude (wrapped around +-180°)
    lon(ii)=wrapToPi(lon(ii));
end


lat=delta;                  %latitude
% put a nan in the lon and lat vectors in case of change from + to - 180°
%to ensure the plot is smooth.

jj=1;
while jj<=length(lon)-1

    if lon(jj)>0 && lon(jj+1)<=0
        lon= [lon(1:jj);pi;nan;-pi;lon(jj+1:end)];
        lat= [lat(1:jj);lat(jj+1);nan;lat(jj);lat(jj+1:end)];
        jj=jj+3;
    end
    jj=jj+1;
    
end

% omega_m=omega_m*pi/(180*3600);
% lon=alpha-th_G0-omega_m*(tspan'-t0)+lon_0;
% lon=wrapToPi(lon);

end