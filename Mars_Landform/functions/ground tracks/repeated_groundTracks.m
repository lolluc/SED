function a = repeated_groundTracks(omegaE, k, m, mu)
%repeated_groundTracks computes the modified semi-major axis of the orbit
%in order to have periodically repeated ground tracks
%
% PROTOTYPE
% a = repeated_groundTracks(omegaE, k, m)
%
% INPUT:
% omegaE [1] Planet's rotation velocity                [deg/h]
% k      [1] Number of revolutions of the satellite    [-]
% m      [1] Number of rotations of the planet         [-]
% mu     [1] Gravitational parameter                   [km^3/s^2]
%
% OUTPUT:
% a      [1] Modified semi-major axis                  [km]
%
% Gives the modified semi-major axis for repeating ground tracks
omegaE = omegaE*pi/(180*3600);

% Computation of the mean anomaly:

n = omegaE*k/m;

% Resulting a from mean anomaly:

a = (mu/n^2)^(1/3);
end

