function [alpha, delta, lon, lat] = groundTrack(kep_el, lon_0, tspan, mu, omegaE, t0, th_G0)
% groundTrack computes the ground track of an orbit
%
% PROTOTYPE
% [alpha, delta, lon, lat] = groundTrack(kep_el, lon_0, tspan, mu, omegaE, t0, th_g0)
%
% INPUT:
% kep_el   [1x6] State of the orbit at initial time (keplerian elements)            [km, -, rad]
% lon_0    [1]   Longitude of Greenwich meridian at initial time                    [rad]
% tspan    [:,1] Vector of times at which the ground track will be computed         [s]
% mu       [1]   Gravitational parameter                                            [km^3/s^2]
% omegaE   [1]   Planet's rotation velocity                                         [deg/h]
% t0       [1]   Reference initial time                                             [s]
% th_G0    [1]   Initial position of the rotating planet with respect to x axis     [rad]
%
% OUTPUT:
% alpha    [:,1] Right ascension in Earth Centered Equatorial Inertial frame        [rad]
% delta    [:,1] Declination in Earth Centered Equatorial Inertial frame            [rad]
% lon      [:,1] Longitude with respect to rotating Earth (0 at Greenwich meridian) [rad]
% lat      [:,1] Latitude with respect to rotating Earth                            [rad]
%
% Gives the ground track of an orbit

[r0,v0]=kep2car(kep_el(1),kep_el(2),kep_el(3),kep_el(4),kep_el(5),kep_el(6),mu);

% Initial conditions
y0 = [r0'; v0'];

% Set options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);

% Perform the integration
[~, Y1] = ode113(@(t, y) ode_orbit(t,y,mu), tspan, y0, options);

% Find state vectors
r_orbit = [Y1(:,1),Y1(:,2),Y1(:,3)];

x=r_orbit(:,1);
y=r_orbit(:,2);
z=r_orbit(:,3);
R = sqrt(x.^2+y.^2+z.^2);
delta=asin(z./R);
lat=delta;
lat=wrapToPi(lat);

alpha=zeros(length(tspan),1);
for j=1:length(tspan)
    if y(j)/R(j)>0
        alpha(j)=acos(x(j)/R(j)/cos(delta(j)));
    else
        alpha(j)=2*pi-acos(x(j)/R(j)/cos(delta(j)));
    end
end

omegaE=omegaE*pi/(180*3600);
lon=alpha-th_G0-omegaE*(tspan'-t0)+lon_0;
lon=wrapToPi(lon);

end

