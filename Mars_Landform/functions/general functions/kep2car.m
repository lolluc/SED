function [r,v] = kep2car(a, e, i, OM, om, th, mu)

% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates
%
% PROTOTYPE:
%   [r, v] = kep2car(a, e, i, OM, om, th, mu)
%
% DESCRIPTION:
%   Conversion from Keplerian elements to Cartesian coordinates.
%
% INPUT:
%   a     [1x1]   Semi-major axis             [km]
%   e     [1x1]   Eccentricity                [-]
%   i     [1x1]   Inclination                 [rad]
%   OM    [1x1]   RAAN                        [rad]
%   om    [1x1]   Pericentre anomaly          [rad]
%   th    [1x1]   True anomaly                [rad]
%   mu    [1x1]   Gravitational parameter     [km^3/s^2]
%
% OUTPUT:
%   r     [1x3]   Position vector             [km]
%   v     [1x3]   Velocity vector             [km/s]

%% define p, r
p     = a*(1-e^2);                              % Semi-latus rectum [km]
r_abs = p/(1+e*cos(th));                        % Absolute value of the position vector [km]

%% define state vector about the perifocal sistem
r_pf  = r_abs*[cos(th), sin(th), 0];            % Position vector (perifocal coordinate system) [km]
v_pf  = sqrt(mu/p)*[-sin(th), e+cos(th), 0];    % Velocity vector (perifocal coordinate system) [km/s]

%% rotate the state vector in the equatorial geocentric sistem

R3_OM = [cos(OM) sin(OM) 0
         -sin(OM) cos(OM) 0
         0 0 1];                                % Rotation matrix for OM around k axis [-]
R1_i  = [1 0 0
         0 cos(i) sin(i)
         0 -sin(i) cos(i)];                     % Rotation matrix for i around i' axis [-]
R3_om = [cos(om) sin(om) 0
         -sin(om) cos(om) 0
         0 0 1];                                % Rotation matrix for om around k'' axis [-]

T_pf_to_eci = R3_OM'*R1_i'*R3_om';              % Global rotation matrix [-]

%% results

r = (T_pf_to_eci*r_pf')';                       % Position vector (Geocentric equatorial coordinates) [km]

v = (T_pf_to_eci*v_pf')';                       % Velocity vector (Geocentric equatorial coordinates) [km/s]

end