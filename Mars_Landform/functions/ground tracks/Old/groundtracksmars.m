function [lon,lat,delta,alpha] =...
    groundtracksmars(a,e,i,OM,om,th_0,t_fin,perturbation,th_gw0,plot_enable)
%
% plots ground-tracks for a specified orbit
%
%
% INPUTS: 
%    a [1x1]       -semi-major axis [km]
%    e [1x1]       -eccentricity [-]
%    i [1x1]       -inclination [rad]
%    OM [1x1]      -RAAN [rad]
%    om [1x1]      -perigee anomaly [rad]
%    th_0 [1x1]    -real anomaly of initial point [rad]
%    t_fin [1x1]   -final time of the propagation [s]
%    perturbation  -string setting of the included perturbation model:
%                          -  'unperturbed'     -use unperturbed 2bp equations.
%                          -  'j2'              -use j2 perturbation model.
%    th_gw0 [1x1]  -angle between greenwich meridian for the earth and gamma direction at
%                   intial time instant. [rad]
%    plot_enable   -string setting to choose the wanted plots:
%                          - 'all'  -plot 3D orbit and ground-tracks.
%                          - 'g-t'  -plot ground-tracks only.
%                          - 'off'  -disables all plots.
%    
%
% OUTPUTS:
%    lon [1xN]     -longitude position history. [rad]
%    lat [1xN]     -latitude position history. [rad]
%    delta [1xN]   -declination position history. [rad]
%    alpha [1xN]   -right ascension position history. [rad]   
%
% REQUIRED FUNCTIONS: astroConstants, solarsystem3D, kep2car
%
%


mu= astroConstants(14); %mars's mu
om_mars=deg2rad((360/24.6229)/3600); % mars's rotational speed


%inital state vector calculated from keplerian elements
[r0,v0]=kep2car(a ,e ,i ,OM,om,th_0);
y0 = [r0;v0];
N_points = 100000;    % Number of points computed
tspan = linspace(0,t_fin,N_points);
R_m = astroConstants(24);
j2 = 0.00196045;
%integration
switch perturbation
    
    case 'unperturbed'
options= odeset('reltol',1e-13,'AbsTol',1e-14);
[t,u]= ode113(@(t, y) ode_orbit(y,mu), tspan, y0, options);
r= u(1:end,1:3);
    case 'j2'
options= odeset('reltol',1e-13,'AbsTol',1e-14);
[t,u]= ode113(@(t, y) ode_perturbed_orbit(y,mu,R_m,j2), tspan, y0, options);
r= u(1:end,1:3);
end


%pre allocation
delta=zeros(length(t),1);
alpha=delta;
th_gw=delta;
lon=delta;


for ii=1:length(t)
    
    delta(ii)= asin(r(ii,3)/sqrt((r(ii,3))^2+(r(ii,2))^2)+(r(ii,1))^2); %declination
    
    alpha(ii) =atan2(r(ii,2),r(ii,1)); % rigth ascension
    
    th_gw(ii)=th_gw0+om_mars*(t(ii)); % teta of greenwich (time)
    th_gw(ii)=wrapTo2Pi(th_gw(ii));
    
    lon(ii)= alpha(ii)-th_gw(ii); %longitude (wrapped around +-180°)
    lon(ii)=wrapToPi(lon(ii));
end

lat= delta; %latitude

% put a nan in the lon and lat vectors in case of change from + to - 180°
%to ensure the plot is smooth.

ii=1;
while ii<=length(lon)-1

    if lon(ii)<0 && lon(ii+1)>=0
        lon= [lon(1:ii);-pi;nan;pi;lon(ii+1:end)];
        lat= [lat(1:ii);lat(ii+1);nan;lat(ii);lat(ii+1:end)];
        ii=ii+3;
    end
    ii=ii+1;
    
end

% plotting

switch plot_enable
    case 'all' %plot both orbit and ground tracks 
        % plot ground tracks
        figure
        I =imread('MarsTexture.jpg');
        image([-180,180],[90 , -90],I)
        hold on
        plot(rad2deg(lon),rad2deg(lat),'g','linewidth',2,'handlevisibility','off')
        plot(rad2deg(lon(1)),rad2deg(lat(1)),'yo','markersize',5,'linewidth',4)
        plot(rad2deg(lon(end)),rad2deg(lat(end)),'rx','markersize',13,'linewidth',4)
        grid on
        xlabel('Longitude [°]')
        ylabel('Latitude [°]')
        leg=legend('Inital point','Final point');
        set(leg, 'TextColor','w', 'Color','none', 'EdgeColor','w')
        ax = gca;
        ax.YDir = 'normal';
        ax.GridColor = 'white';
        
        % plot 3D orbit
        solarsystem3D(4,'black')
        plot3(r(1:end,1),r(1:end,2),r(1:end,3),'color',[23 115 238]/255,'linewidth',2)
        plot3(r(1,1),r(1,2),r(1,3),'rx','linewidth',4)
        plot3(r(end,1),r(end,2),r(end,3),'yx','linewidth',4)
        
    case 'g-t' %plot only ground tracks
         figure
        I =imread('MarsTexture.jpg');
        image([-180,180],[90 , -90],I)
        hold on
        plot(rad2deg(lon),rad2deg(lat),'g','linewidth',2,'handlevisibility','off')
        plot(rad2deg(lon(1)),rad2deg(lat(1)),'yo','markersize',5,'linewidth',4)
        plot(rad2deg(lon(end)),rad2deg(lat(end)),'rx','markersize',13,'linewidth',4)
        grid on
        xlabel('Longitude [°]')
        ylabel('Latitude [°]')
        leg=legend('Inital point','Final point');
        set(leg, 'TextColor','w', 'Color','none', 'EdgeColor','w')
        ax = gca;
        ax.YDir = 'normal';
        ax.GridColor = 'white';
end
end


