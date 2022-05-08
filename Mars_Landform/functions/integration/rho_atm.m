function [rho]= rho_atm(h)
% RHO_ATM Exponential atmospheric model
%
% INPUTS:
%   h[1x1]    - height from surface [km]
%
% OUTPUTS:
%   rho[1x1]  - air density [kg/m^3]
% 
%SOURCE: wertz

if h<=0
    rho0=0 ;H=1;     h0=0;
elseif h>=1000
    h0=1000; rho0=3.019e-15; H=268.00 ;
elseif h<25 
    h0=0 ; rho0=1.225; H=7.249;
elseif h<30
    h0=25 ; rho0=3.899e-2; H=6.349;
elseif h<40 
    h0=30 ; rho0=1.774e-2; H=6.682;
elseif h<50 
    h0=40 ; rho0=3.972e-3; H=7.554;    
elseif h<60 
    h0=50 ; rho0=1.057e-3; H=8.382;    
elseif h<70 
    h0=60 ; rho0=3.206e-4; H=7.714;    
elseif h<80 
    h0=70 ; rho0=8.770e-5; H=6.549;    
elseif h<90 
    h0=80 ; rho0=1.905e-5; H=5.799;    
elseif h<100 
    h0=90 ; rho0=3.396e-6; H=5.382;    
elseif h<110 
    h0=100 ; rho0=5.297e-7; H=5.877;    
elseif h<120 
    h0=110 ; rho0=9.661e-8; H=7.263;    
elseif h<130 
    h0=120 ; rho0=2.438e-8; H=9.473;    
elseif h<140
    h0=130; rho0=8.484e-9; H=12.636 ;  
elseif h<150 
    h0=140; rho0=3.845e-9; H=16.149 ; 
elseif h<180
    h0=150; rho0=2.070e-9; H=22.523 ;  
elseif h<200 
    h0=180; rho0=5.464e-10; H=29.740 ; 
elseif h<250
    h0=200; rho0=2.789e-10; H=37.105 ;  
elseif h<300 
    h0=250; rho0=7.248e-11; H=45.546 ; 
elseif h<350
    h0=300; rho0=2.418e-11; H=53.628 ;  
elseif h<400 
    h0=350; rho0=9.158e-12; H=53.298 ; 
elseif h<450
    h0=400; rho0=3.725e-12; H=58.515 ;  
elseif h<500 
    h0=450; rho0=1.585e-12; H=60.828 ; 
elseif h<600
    h0=500; rho0=6.967e-13; H=63.822 ;  
elseif h<700 
    h0=600; rho0=1.454e-13; H=71.835 ; 
elseif h<800
    h0=700; rho0=3.614e-14; H=88.667 ;  
elseif h<900 
    h0=800; rho0=1.170e-14; H=124.64 ;     
elseif h<1000
    h0=900; rho0=5.245e-15; H=181.05 ;    
end


rho=rho0*exp(-(h-h0)/H);

return