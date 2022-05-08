function [deltat] = fligthtime(a,e,teta1,teta2,mu)
% DESCRIPTION:
% Calculate time between known real anomalies on a specified orbit.
%
% INPUT:
% a [1x1] Semi-major axis [km]
% e [1x1] Eccentricity [-]
% teta1 [1x1] True anomaly point 1 [rad]
% teta2 [1x1] True anomaly point 2 [rad]
% mu [1x1] Gravitational parameter [km^3/s^2]
%           (default value of Earth, if you have a different focal centre
%           please specify the value)
%
% OUTPUT:
% deltat [1x1] time period between the two points [s]
if nargin<5
    mu = 398600.433 ;
end


if a>0
    
    
    n=sqrt(mu/a^3); T= 2*pi/n;
    
    if teta1>2*pi
        error('theta 1 must be less than 2pi')
    end
    k=0;
    while teta2>2*pi
        teta2=teta2-2*pi;
        k=k+1;
    end
    %% define times at points 1 & 2
    
    t1=(1/n)*(1/(1-e^2)^(3/2))*(2*atan(sqrt((1-e)/(1+e))*tan(teta1/2))-...
        e*sqrt(1-e^2)*sin(teta1)/(1+e*cos(teta1)));
    
    t2=(1/n)*(1/(1-e^2)^(3/2))*(2*atan(sqrt((1-e)/(1+e))*tan(teta2/2))-...
        e*sqrt(1-e^2)*sin(teta2)/(1+e*cos(teta2)));
    
    %% considerationts about the signum of the results
    
    if t1<0
        t1=T+t1;
    end
    if t2<0
        t2=T+t2;
    end
    
    %% calculate time period
    
    deltat= t2-t1;
    if deltat<0
        deltat=T+deltat;
    end
    
    deltat=deltat+k*T;
    
else
    theta_inf=acos(-1/e);
    if teta2>theta_inf
        error('theta>theta_inf')
    end
    E=2*atanh(sqrt((e-1)/(e+1))*tan(teta1/2));
    M=e*sinh(E)-E;
    t1=M*sqrt(-a^3/mu);
    E=2*atanh(sqrt((e-1)/(e+1))*tan(teta2/2));
    M=e*sinh(E)-E;
    t2=M*sqrt(-a^3/mu);
    deltat=abs(t1-t2);
    
    
end
    
    return