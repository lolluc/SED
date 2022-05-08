%% Lab 2
close all
clear
clc

%% Data

mu = astroConstants(14);
a  = astroConstants(24) + 100;
e  = 0;
i  = 80*pi/180;
OM = 0*pi/180;
om = 0*pi/180;
th0 = 0*pi/180;
kep_el = [a e i OM om th0];
lon_0 = 0;      % Initial longitude
T = 2*pi*sqrt(a^3/mu);
N = 50;         % Number of orbits
N_points = 100000;    % Number of points computed
tspan = linspace(0,N*T,N_points);
omegaE = 14.62;  % Planet's angular velocity
t0 = 0;      % Initial time
th_G0 = 0;   % Theta_G at initial time
% Repeating parameters
k = 1;
m = 1;
% Mars parameters
Re = astroConstants(24);
j2 = 0.00196045;

%% Exercise 1

% Computation of the ground track
[~, ~, lon, lat] = groundTrack(kep_el, lon_0, tspan, mu, omegaE, t0, th_G0);

% Plot of the ground track
for ii=2:length(lon)
    if abs(lon(ii))>0.5
        if lon(ii-1)*lon(ii)<0
            lon(ii-1)=NaN;
        end
    end
end

lon=lon*180/pi;
lat=lat*180/pi;
figure
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon,lat,'g','linewidth',1.2);
plot(lon(1),lat(1),'go','linewidth',2)
plot(lon(length(lon)),lat(length(lat)),'gs','linewidth',2)
xlim([-180 180])
ylim([-90 90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('Ground track','Start','End','Location','best','Orientation','horizontal')
title('Ground track plot')

%% Exercise 2

% Plot of the original ground track

figure
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon,lat,'g','linewidth',1.2);
xlim([-180 180])
ylim([-90 90])

% Modification of the semi-major axis

a_new = repeated_groundTracks(omegaE, k, m, mu);
kep_el(1) = a_new;
T_new = 2*pi*sqrt(a_new^3/mu);
tspan_new = linspace(0,N*T_new,N_points);

% Computation of the modified ground track

[~, ~, lon_rep, lat_rep] = groundTrack(kep_el, lon_0, tspan_new, mu, omegaE, t0, th_G0);

% Plot of the modified ground track for confrontation

for ii=2:length(lon_rep)
    if abs(lon_rep(ii))>0.5
        if lon_rep(ii-1)*lon_rep(ii)<0
            lon_rep(ii-1)=NaN;
        end
    end
end

lon_rep=lon_rep*180/pi;
lat_rep=lat_rep*180/pi;
plot(lon_rep,lat_rep,'r','linewidth',1.2);
plot(lon(1),lat(1),'go','linewidth',2)
plot(lon(length(lon)),lat(length(lat)),'gs','linewidth',2)
plot(lon_rep(1),lat_rep(1),'ro','linewidth',2)
plot(lon_rep(length(lon_rep)),lat_rep(length(lat_rep)),'rs','linewidth',2)
xlim([-180 180])
ylim([-90 90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('Original','Repeating','Start','End','Start','End','Location','best','Orientation','horizontal')
title('Repeating ground track confrontation')

%% Exercise 3

% Modification of the semi-major axis

a_perturbed = j2_groundTrack(omegaE, k, m, mu, j2, Re, a, kep_el(2), kep_el(3));
kep_el(1) = a_perturbed;
T_perturbed = 2*pi*sqrt(a_perturbed^3/mu);
tspan_perturbed = linspace(0,N*T_perturbed,N_points);

% Computation of the perturbed ground track
[~, ~, lon_per, lat_per] = perturbed_groundTrack(kep_el, lon_0, tspan_perturbed, mu, omegaE, t0, th_G0, Re, j2);

% Plot of the original ground track
kep_el(1) = a;
[~, ~, lon1, lat1] = perturbed_groundTrack(kep_el, lon_0, tspan, mu, omegaE, t0, th_G0, Re, j2);
for ii=2:length(lon1)
    if abs(lon1(ii))>0.5
        if lon1(ii-1)*lon1(ii)<0
            lon1(ii-1)=NaN;
        end
    end
end

lon1=lon1*180/pi;
lat1=lat1*180/pi;
figure
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon1,lat1,'g','linewidth',1.2);
xlim([-180 180])
ylim([-90 90])

% Plot of the repeating unperturbed ground track
kep_el(1) = a_new;
[~, ~, lon2, lat2] = perturbed_groundTrack(kep_el, lon_0, tspan_new, mu, omegaE, t0, th_G0, Re, j2);
for ii=2:length(lon2)
    if abs(lon2(ii))>0.5
        if lon2(ii-1)*lon2(ii)<0
            lon2(ii-1)=NaN;
        end
    end
end

lon2=lon2*180/pi;
lat2=lat2*180/pi;
plot(lon2,lat2,'r','linewidth',1.2);

% Plot of the repeating ground track for secular j2

for ii=2:length(lon_per)
    if abs(lon_per(ii))>0.5
        if lon_per(ii-1)*lon_per(ii)<0
            lon_per(ii-1)=NaN;
        end
    end
end

lon_per=lon_per*180/pi;
lat_per=lat_per*180/pi;
plot(lon_per,lat_per,'y','linewidth',1.2);
plot(lon1(1),lat1(1),'go','linewidth',2)
plot(lon1(length(lon1)),lat1(length(lat1)),'gs','linewidth',2)
plot(lon2(1),lat2(1),'ro','linewidth',2)
plot(lon2(length(lon2)),lat2(length(lat2)),'rs','linewidth',2)
plot(lon_per(1),lat_per(1),'yo','linewidth',2)
plot(lon_per(length(lon_per)),lat_per(length(lat_per)),'ys','linewidth',2)
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('Original','Repeating GT for unperturbed','Repeating GT for secular j2', ...
       'Start','End','Start','End','Start','End','Location','best')
title('Repeating ground track confrontation')
