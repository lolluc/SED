function a = j2_groundTrack(omegaE, k, m, mu, j2, Re, a0, e, i, tol)
%repeated_groundTracks computes the modified semi-major axis of the orbit
%in presence of j2 effect
%
% PROTOTYPE
% a = j2_groundTrack(omegaE, k, m, mu, j2, Re, a0, e, i, tol)
%
% INPUT:
% omegaE [1] Planet's rotation velocity                [deg/h]
% k      [1] Number of revolutions of the satellite    [-]
% m      [1] Number of rotations of the planet         [-]
% mu     [1] Gravitational parameter                   [km^3/s^2]
% j2     [1] Perturbation parameter                    [-]
% Re     [1] Radius of the planet                      [km]
% a0     [1] Semi-major axis non perturbed             [km]
% e      [1] Eccentricity                              [-]
% i      [1] Inclination                               [rad]
% tol    [1] Tolerance of the intergation              [-]
%
% OUTPUT:
% a      [1] Modified semi-major axis                  [km]
%
% Gives the modified semi-major axis for perturbed orbits
omegaE = omegaE*pi/(180*3600);

d_OM = @(a) -(3/2*(sqrt(mu)*j2*Re^2)/((1-e^2)^2*a^(7/2)))*cos(i);

d_om = @(a) -(3/2*(sqrt(mu)*j2*Re^2)/((1-e^2)^2*a^(7/2)))*(5/2*(sin(i))^2-2);

d_M0 = @(a) 3/2*(sqrt(mu)*j2*Re^2)/((1-e^2)^(3/2)*a^(7/2))*(1-3/2*(sin(i))^2);

if nargin<10
    tol = 1e-6;
end
options = optimset('TolX',tol,'Display','off');

a = fzero(@(a) (m/k)-(omegaE-d_OM(a))/(sqrt(mu/a^3)+d_om(a)+d_M0(a)),a0,options);

end

