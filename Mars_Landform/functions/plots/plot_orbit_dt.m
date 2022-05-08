function [Rxo, Ryo, Rzo] = plot_orbit_dt (a, e, i, OMEGA, omega, th0, thf, dth, mu)

%Plot orbita 3d

% Input:
% a         [1x1]   semiasse maggiore           [km]
% e         [1x1]   eccentricità                [-]
% i         [1x1]   inclinazione                [rad]
% OMEGA     [1x1]   RAAN                        [rad]
% omega     [1x1]   anomalia del pericentro     [rad]
% th0       [1x1]   anomalia vera iniziale      [rad]
% thf       [1x1]   anomalia vera finale        [rad]
% dth       [1x1]   passo anomalia vera         [rad]
% mu        [1x1]   parametro gravitazionale    [km^3/s^2]


if th0 < thf 
    theta = [th0 : dth : thf];                 %definisco divisione in punti orbita
else 
    theta = [th0 : dth : thf + 2*pi];
end 

p = (1-e^2) * a ;                               %Definisco p per l'orbita 
erre = @(theta) p ./ ( 1 + e*cos(theta) );      %Scrivo l'equazione dell'orbita (ellisse) rispetto al sistema di riferimento perifocale
R= [erre(theta)];                              %Definisco un vettore che descrive l'orbita calcolando l'equazione nei punti alpha

x= [R.*cos(theta)];                             %Scompongo R nelle componenti rispetto al sistema perifocale
y=[R.*sin(theta)];
z= [R.*zeros];

R_pf= [x; y; z];

% Trasformazione da sistema perifocale a sistema geocentrico orbita 
[Rot_ge_pf,Rot_pf_ge] = rot (OMEGA,omega,i);

R_ge= Rot_pf_ge*R_pf;

% Plot orbita 
Rxo=R_ge(1,:);
Ryo=R_ge(2,:);
Rzo=R_ge(3,:);


end