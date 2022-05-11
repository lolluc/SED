%% CODE FOR TELECOMMUNICATION SYSTEMS DESIGN - ARCAS proejct, comp_telecom.m 

%  DESCRIPTION:
%  This is the script for the design of a telecommunication system. 
%  In particular, this code is capable to check the consistency of a set of
%  antennas for the communication between a lander, on the Callisto Moon (of
%  Jupiter), and the antennas of the Europa Clipper NASA space mission
%  during its flybys on Callisto.
%  -----------------------------------------------------------------------
%  
%  AUTHOR: Davide Gravina (davide.gravina@mail.polimi.it)
%  DATE: 30 June 2021 (last version)
% 
%  ----------------------------------------------------------------------- 
%  SECTIONS:
%  1. Computation of Power Density (2D model) of lander during the orbiter flybys.
%  2. Computation of Power Density (3D model) of lander during the orbiter flybys.
%  3. Computation of Losses.
%  4. Computation of Received Power.
%  5. Computation of Noise.
%  6. Computation of Link Budget.
%  7. Computation of effective Signal to Noise Ratio (SNR).
%  8. Computation of Data rate.
%  9. Selection of Proper Bandwidth.
%  10.Trials...

clear all
close all

%% COMPUTATION OF POWER DENSITY (2D MODEL):
clear all
close all

% Data:
R_cal=2410.3; %[Km] radius of Callisto Moon
h_min=25; %[Km] minimum altitude of flybys
h_max=2525; %[Km] maximum altitude of flybys
h=[h_min:100:h_max]; %[km] vector of all altitude
el=[0:1:180]; %[deg] elevation angle
lat=[-90:1:90]; %[deg] latitude of flyby;

% NOTE: Europa clipper have more antennas in order to achieve a complete
% pattern of sending and reciving area along all the y-z plane (x=
% direction of satellite movement/orbit).

% 2D triangle:
a=R_cal+h;
b=R_cal;
alpha=el+90;
theta=lat;
%(find c, and limiting parameters)

% % Computation of max latitude with elevation angle=0:
% %(Assume to be placed at equator);
% c1=sqrt(a.^2-b.^2); %[km] distance of flyby from lander
% lat_max=rad2deg(asin(c1./a)); %[deg] max latitude with elevation angle=0

% Computation of max latitude, GENERAL WAY: (CHECK=OK)
%(Assume to be placed at equator, valid for any point considering different initial latitude)

%Choose a fixed elevation angle:
el1_min=10; %[deg] minimum elevation angle, considered.
alpha1_min=90+el1_min; %[deg] 

%Values of triangle for minimum elevation angle (at each possible altitude)
beta2=rad2deg(asin((b./a)*sin(deg2rad(alpha1_min)))); %[deg] angle between lander-sat (link) and sat-ground
c2= a.*cos(deg2rad(beta2)) + b.*cos(deg2rad(alpha1_min)); %[km] distance of link lander-sat
lat_max2= rad2deg(asin((c2./a).*sin(deg2rad(alpha1_min)))); %[deg] maximum latitude for fixed elevation angle at each altitude.

%Results: 2D model
% % % figure()
% % % [h_x, Lat_y1, Distance_y2]= plotyy(h,lat_max2,h,c2);
% % % title('Values with minimum Elevation-Angle=5°');
% % % xlabel('Altitude of Flybys [km]');
% % % ylabel(h_x(1),'Maximum Latitude [deg]');
% % % ylabel(h_x(2),'Distance of signal [km]');
% % % grid on

% figure()
% 
% subplot(1,2,1)
% plot(h,lat_max2,'r','LineWidth',1.5);
% title('2D : Values with minimum Elevation-Angle=5°');
% xlabel('Altitude of Flybys [km]');
% ylabel('Maximum Latitude [deg]');
% grid on
% 
% subplot(1,2,2)
% plot(h,c2,'b','LineWidth',1.5);
% title('2D : Values with minimum Elevation-Angle=5°');
% xlabel('Altitude of Flybys [km]');
% ylabel('Maximum travel of signal [km]');
% grid on
% 
% 

%% COMPUTATION OF POWER DENSITY (3D MODEL):
%Semplification, instead using the spherical geometry it will be use
%planare geometry for computation of distance of signal in all the cases,
%with an reliable approximation (Usando la dioganale di un parallelepipedo
%retto dovrebbe essere lo stesso identico valore).

%Assuming as maximum value, due to a minimum elevation angle, the value of
%max_lat2, also for longitude, it is computed the new maximum distance of
%signal for a 3D case (so with a satellati that it's not on the plane of
%lander.
% initial_lat=0;%[deg]
% initial_long=0;%[deg]
max_Delta_latitude= lat_max2; %[deg] +- this value from the initial lat of lander
max_Delta_longitude= lat_max2; %[deg] +- this value from the initial lat of lander
%THE USE OF lat_max2 FOR BOTH THE EQUATION IS CORRECT, BECAUSE THERE ARE THE SAME.

%Planar Approsimation:
x_S=2*a.*sin(deg2rad(max_Delta_longitude./2)); %[km]
y_S=2*a.*sin(deg2rad(max_Delta_latitude./2)); %[km]

new_c2=sqrt((x_S.^2) + (y_S.^2) + (h.^2)); %[km]

%There is one APPROXIMATION: the maximum distance is computed in the
%pericentre of flyby, but after (and before) that point the s/c is faraway
%respect the pericentre position. _is OK 
% (not so different a part for the limit case of h=2500 km)
% {For Example: h_p=25 km -> at True anomaly=52° := r=748 km}

%Result: 3D model
figure()

subplot(1,3,1)
plot(h,lat_max2,'r','LineWidth',1.5);
title('3D : Values with minimum Elevation-Angle=10°');
xlabel('Altitude of Flybys [km]');
ylabel('Maximum Latitude [deg]');
grid on

subplot(1,3,2)
plot(h,lat_max2,'r','LineWidth',1.5);
title('3D : Values with minimum Elevation-Angle=10°');
xlabel('Altitude of Flybys [km]');
ylabel('Maximum Longitude [deg]');
grid on

subplot(1,3,3)
plot(h,new_c2,'b','LineWidth',1.5);
title('3D : Values with minimum Elevation-Angle=10°');
xlabel('Altitude of Flybys [km]');
ylabel('Maximum travel of signal [km]');
grid on
%% Comparison 2D/3D
% figure()
% subplot(1,2,1)
% plot(h,c2,'b','LineWidth',1.5);
% title('2D : Values with minimum Elevation-Angle=5°');
% xlabel('Altitude of Flybys [km]');
% ylabel('Maximum travel of signal [km]');
% grid on
% subplot(1,2,2)
% plot(h,new_c2,'b','LineWidth',1.5);
% title('3D : Values with minimum Elevation-Angle=5°');
% xlabel('Altitude of Flybys [km]');
% ylabel('Maximum travel of signal [km]');
% grid on



%% COMPUTATION OF LOSSES
c=299792.456; %[km/s] speed of light
r= linspace(new_c2(1),new_c2(end),1000); %[km] maximum distances of signal 
f=8.15; %[GHz] frequencies of signal - X-Band
% X-band typical values: UPLINK= 7.9/8.4 GHz ; DOWNLINK 7.25/7.75 GHz
wl=(c*1e3)/(f*1e9); %[m]=[m/s]/[Hz] Wavelength of (f) frequency in Free Space (FS)


%Free space losses:
L_s= 20.*log10(r.*f) + 92.45; %[dB] (from wikipeda)_is OK
%WHERE: r = distance in [km]; f= frequency in [GHz]; (+92.45 = Constant = 20*log(4*pi/c) with [particular units of measure])


%Misalignement losses: 
%ASSUME: Perfect direction of link and simplified computation
L_mis=100; %[%] PErcentage of correct direction of link(1%<L_mis<100%)
L_mis=10*log10(L_mis/100); %[dB] ASSUME a perfect directional link


%Atmospheric Losses: RESULT: NO LOSS DUE TO ATTENUATION !!!

%Evanescent wave problem:
N_max=2e10 ; %[e-/m^3] electron density peak of Callisto (from 1e10 to 2e10)
f_plasma=9*sqrt(N_max); %[Hz] plasma frequency = limit for evanescent wave
% % f_plasma_max=1.2728e+06 [Hz] =1.27 MHz - we are quite above this limit.

%Permittivity (epsi) of Callisto Ionosphere:
e_charge= 1.6e-19; %[Coulomb = C] electron charge
e_mass= 9.10e-31; %[kg] electron mass
N_max=2e10; %[e-/m^3] electron density peak of Callisto (from 1e10 to 2e10)
af=2*pi*(f*1e9); %[rad/s] ANGULAR FREQUENCY of signal frequency

epsi_0= 8.85418e-12; %[F/m] permittivity of free space (or air)
epsi_rel= 1 - ((N_max*(e_charge^2))/((af^2)*epsi_0*e_mass)); % relativi permittivity due to ionosphere of Callisto
%epsi_rel= 1, so FREE SPACE = ininfluent ionosphere_is OK

epsi = epsi_0*epsi_rel; %[F/m] equal to Permittivity of Free Space (o air) 

%So, also Permeability (perm) of Callisto is equal to Permability of FreeSpace (o Air)
perm_0=1.256637e-6; %[H/m]
perm_rel=1; % equal to Free Space
perm = perm_0*perm_rel; %[H/m]



%% Power Receiver:

%REMEMBER: conversion to dB: {https://www.researchgate.net/post/What-is-the-formula-for-converting-decibels-into-amplitude-magnitude#:~:text=a%20same%20scale.-,For%20converting%20the%20ratio%20of%20two%20power%20values%20in%20dB,ans(dB)%2F20).}
% For converting the ratio of two power values in dB, we use ans(dB) = 10*log10(ratio) and ratio=10^(ans(dB)/10).
% For converting the ratio of two voltage or current values in dB, we use ans(dB)=20*log10(ratio) and ratio=10^(ans(dB)/20). 

P_tx=[1:1:4]; %[W] power transmitted from Trasmitter=Lander (I can choose in this range)
P_tx_dB= 10*log10(P_tx); %[dB]-(ratio computed respect 1W) power transmitted from Trasmitter=Lander

%TRASMITTER 1 - Planar array antenna: 
% % https://www.endurosat.com/cubesat-store/cubesat-antennas/x-band-4x4-patch-array/
% (Diameter chosed arbitrarly= a circle of Diameter of D, that rapresent a square of 10cm (Planar Array))
% D_tx=0.12;%[m] Diameter of antenna (tx: lander)
% Area_tx= 0.01; %[m^2] Area of Trasmitter = Planar Array antenna (SQUARE= 10x10 cm)
% effic_tx=0.60; % (ASSUMPTION) efficiency of antenna (real efficiency for antenna (form 0.5 to 0.6)_ is OK
% effic_dB_tx=10*log10(effic_tx); %[dB] efficiency of antenna expressed in deciBel
% % Peak Gain [dBi]: (where wavelength(wl)=0.0368[m], equal to f=8.15GHz)
% pG_tx =10*log10(Area_tx/(wl^2))+8; %[dBi]= deciBel isotropic - Peak Gain of Planar Array antenna
% % Directivity [dB]: equal to Peak Gain, but with another unit of measure
% Dm_tx =pG_tx-2.15; %[dB]= deciBel Directivity of Planar Array antenna
% % Gain trasmitter [dB]:
% G_tx_dB=effic_dB_tx+Dm_tx; %[dB]
% % G_tx_dB= 16.15 -2.15; %[dB] fixed value from datasheet 
% % % HPBW(Half Power Beam Width) : (Parabola equation)
% % D_tx= 2*sqrt(Area_tx/pi); %[m] Equivalent area of a circle (that approximate our square)
% % % theta_3dB_tx=65*wl/D_tx; %[deg] computed approximated value of HPBW(Half Power Beam Width)
% % theta_3dB_tx= 18; %[deg] fixed value from datasheet 


% % Parabolic Antenna: (Two possible formulas) TRASMITTER 2
Diam_tx= 0.7; %[m] Diameter of parabolic Trasmitter = Parabola antenna (equivalent diameter= 12 cm)
effic_tx=0.60; % (ASSUMPTION) efficiency of antenna (real efficiency for antenna (form 0.5 to 0.6)_ is OK
effic_dB_tx=10*log10(effic_tx); %[dB] efficiency of antenna expressed in deciBel
% Peak Gain [dBi]: TWO POSSIBLE FORMULA - (where wavelength(wl)=0.0368[m], equal to f=8.15GHz)
pG_tx1 = 10*log10(((pi*Diam_tx)/wl)^2); %[dB] Peak Gain of Parabola antenna
pG_tx2 = 20*log10(f)+20*log10(Diam_tx)+17.8; %[dBi]= deciBel isotropic - Peak Gain of Parabola antenna
% Directivity [dB]: equal to Peak Gain, but with another unit of measure
Dm_tx1 =pG_tx1; %[dB]= deciBel Directivity of Parabola antenna - 1st FORMULA
Dm_tx2 =pG_tx2-2.15; %[dB]= deciBel Directivity of Parabola antenna - 2nd FORMULA
% Gain trasmitter [dB]:
G_tx_dB= effic_dB_tx+Dm_tx1; %[dB] Gain Parabolic antenna-1st FORMULA
G_tx2= effic_dB_tx+Dm_tx2; %[dB] gain Parabolic antenna-2nd FORMULA
% (Remember that in dB you sum the two effects, and without dB it's a moltiplication!)
% HPBW(Half Power Beam Width) :
theta_3dB_tx=65*wl/Diam_tx; %[deg] computed value of HPBW(Half Power Beam Width)

% % TRANSMITTER - Omnidirectional antenna (for emergencies)
% G_tx_dB=0; %[dB] typical values 0dBi < G < 3dBi
% theta_3dB_tx= 360; %[deg] it's omnidirectional 

% RECIEVERS- (HGA-High Gain Antenna, MGA-Medium Gain Antenna, LGA- Low Gain Antenna)
%ASSUME: all the values are taken from: https://www.sciencedirect.com/topics/engineering/antenna-gain
% (and from https://didattica-2000.archived.uniroma2.it//ap2/deposito/Day_1.pdf)
% (and from CL#19-0042.pdf =EUROPA CLIPPER)

%[REMEMBER: dB=dBi-2.15]
%LGA- Low Gain Antenna:
G_rx1_dB= 9 -2.15; %[dB] typical values: < 10dBi
theta1_3dB_rx= 180; %[deg] typical values: from Europa Clipper 

%MGA- Low Gain Antenna:
G_rx2_dB= 19 -2.15; %[dB] typical values: 10dBi < G < 20dBi
theta2_3dB_rx= 16; %[deg] typical values: from Europa Clipper 

%MGA- Low Gain Antenna:
G_rx3_dB= 44 -2.15; %[dB] typical values: 20dBi < G < 45dBi
theta3_3dB_rx= 3; %[deg] typical values: from Europa Clipper 


%Radiation Patterns: 
% ASSUME: we are able to have an almost perfect link in Tx and Rx.
f_tx=0.99; %[0<f<1]
f_tx_dB=10*log10(f_tx); %[dB] 
f_rx=0.99; %[0<f<1]
f_rx_dB=10*log10(f_rx); %[dB] 

%Losses of Antenna Devices, for Tx and Rx:
%ASSUME: null losses due the very small device into satellites
L_tx=1; %[0<L<1]
L_tx_dB=10*log10(L_tx); %[dB]
L_rx=1; %[0<L<1]
L_rx_dB=10*log10(L_rx); %[dB]

%MISALLIGNMENT LOSSE:
%Pointing Losses of antenna Tx and Rx: due to Pointing Accuracy
%ASSUME: the accuracy of pointing is very precise theta_accuracy=0.1 [deg]
theta_accuracy_tx= 1.0; %[deg]
% theta_3dB_tx= ...[deg] already computed
L_ptx_dB=-12*((theta_accuracy_tx/theta_3dB_tx)^2);%[dB]
theta_accuracy_rx1=0.01; %[deg] LGA
theta_accuracy_rx2=3; %[deg] MGA
theta_accuracy_rx3=0.1; %[deg] HGA
% theta_3dB_rx= ...[deg] already computed for 3 different possibilities
L_prx1_dB=-12*((theta_accuracy_rx1/theta1_3dB_rx)^2); %[dB] LGA
L_prx2_dB=-12*((theta_accuracy_rx2/theta2_3dB_rx)^2); %[dB] MGA
L_prx3_dB=-12*((theta_accuracy_rx3/theta3_3dB_rx)^2); %[dB] HGA

% Attenuation parameter: 
% Due to the values of ionosphere of Callisto and to the almost (FS)Free
% Space Environment there are NO additional Attenuation effects.
% Maybe there could be some attenuation effects do to the presence of
% Water(H2O) in the atmosphere.
% ASSUME: it's assumed Almost NULL Additional Attenuation
A=0.99; %[0<A<1]
A_dB=10*log10(A); %[dB]

%Already computed Losses:

%Free space losses: [row=1; column=26], where 26 are the (h)altitude_flyby
% L_s= 20.*log10(r.*f) + 92.45; %[dB] 
% WHERE: r = distance in [km]; f= frequency in [GHz]; (+92.45 = Constant = 20*log(4*pi/c) with [particular units of measure])

%Misalignement losses:
% theta_3dB=70*wl/D_tx; %[deg] computed value of HPBW(Half Power Beam Width)
% theta_3dB= 18.5; %[deg] fixed value from datasheet +0.5 deg = error
% L_mis=-12.*((exp(1)/theta_3dB)^2); %[dB] 

for j=1:1:length(P_tx_dB)
%Power Received by Receiver: [dB] Receiver_1 =LGA
P_rx1_dB(j,:)=P_tx_dB(j)+G_tx_dB+f_tx_dB+L_tx_dB+L_ptx_dB + (-L_s+L_mis+A_dB) + G_rx1_dB+f_rx_dB+L_rx_dB+L_prx1_dB;

%Power Received by Receiver: [dB] Receiver_2 =MGA
P_rx2_dB(j,:)=P_tx_dB(j)+G_tx_dB+f_tx_dB+L_tx_dB+L_ptx_dB + (-L_s+L_mis+A_dB) + G_rx2_dB+f_rx_dB+L_rx_dB+L_prx2_dB;

%Power Received by Receiver: [dB] Receiver_3 =HGA
P_rx3_dB(j,:)=P_tx_dB(j)+G_tx_dB+f_tx_dB+L_tx_dB+L_ptx_dB + (-L_s+L_mis+A_dB) + G_rx3_dB+f_rx_dB+L_rx_dB+L_prx3_dB;
end

%Received Power [W]:
P_rx_1=10.^(P_rx1_dB(:,:)./10); %[W] Receiver_1 =LGA
P_rx_2=10.^(P_rx2_dB(:,:)./10); %[W] Receiver_2 =MGA
P_rx_3=10.^(P_rx3_dB(:,:)./10); %[W] Receiver_3 =HGA


%Results: MESSI DOPO (AFTER COMPUTATION OF NOISE)
% figure()
% 
% subplot(1,3,1)
% plot(r,P_rx1_dB(:,:),'LineWidth',1.2)
% title('Power_{Rx} [dB] for LGA Receiver');
% xlabel('Maximum travel of signal [km]');
% ylabel('Recived Power [dB]');
% legend('P_tx=1W','P_tx=2W','P_tx=3W','P_tx=4W','Location','northeast')
% grid on
% 
% subplot(1,3,2)
% plot(r,P_rx2_dB(:,:),'LineWidth',1.2)
% title('Power_{Rx} [dB] for MGA Receiver');
% xlabel('Maximum travel of signal [km]');
% ylabel('Recived Power [dB]');
% legend('P_tx=1W','P_tx=2W','P_tx=3W','P_tx=4W','Location','northeast')
% grid on
% 
% subplot(1,3,3)
% plot(r,P_rx3_dB(:,:),'LineWidth',1.2)
% title('Power_{Rx} [dB] for HGA Receiver');
% xlabel('Maximum travel of signal [km]');
% ylabel('Recived Power [dB]');
% legend('P_tx=1W','P_tx=2W','P_tx=3W','P_tx=4W','Location','northeast')
% grid on


% figure ()
% subplot(1,3,1)
% plot(r,P_rx_1(:,:),'LineWidth',1.2)
% title('Power_{Rx} [W] for LGA Receiver');
% xlabel('Maximum travel of signal [km]');
% ylabel('Recived Power [W]');
% legend('P_tx=1W','P_tx=2W','P_tx=3W','P_tx=4W','Location','northeast')
% grid on
% 
% subplot(1,3,2)
% plot(r,P_rx_2(:,:),'LineWidth',1.2)
% title('Power_{Rx} [W] for MGA Receiver');
% xlabel('Maximum travel of signal [km]');
% ylabel('Recived Power [W]');
% legend('P_tx=1W','P_tx=2W','P_tx=3W','P_tx=4W','Location','northeast')
% grid on
% 
% subplot(1,3,3)
% plot(r,P_rx_3(:,:),'LineWidth',1.2)
% title('Power_{Rx} [W] for HGA Receiver');
% xlabel('Maximum travel of signal [km]');
% ylabel('Recived Power [W]');
% legend('P_tx=1W','P_tx=2W','P_tx=3W','P_tx=4W','Location','northeast')
% grid on


%% NOISE COMPUTATION:

%TYPICAL RANGE OF NOISE FIGURE(NF) FOR A RECIEVER: (NF=2dB/10dB);
%So, the typical equivalent noise temperatures of Reciever (Rx) are: (T_eq= 46K/290K);

%Callisto temperature are: for Lander-Rx(from pdf: "planet and moon properties")
T_min_surface=80; %[K]
T_max_surface=160; %[K]
T_mean_surface=116; %[K]
T_env_rx1= T_mean_surface; %[K] Temperature of environment near Rx antenna - Case 1 Uplink-Lander
% T_env_rx1= T_min_surface; %[K] BEST CASE Temperature of environment near Rx antenna - Case 1 Uplink-Lander

%Space Temperature: for Mothership-Rx
T_env_rx2= 3; %[K] Temperature of environment near Rx antenna - Case 2 Downlink-Lander

%Antenna temperatures: (from Europa clipper datasheet)
%for HGA:
T_min_ant_HGA= 73; %[K]
T_max_ant_HGA= 373; %[K]
T_mean_ant_HGA= 223; %[K]
%for MGA-LGA:
T_min_ant_MLGA= 138; %[K]
T_max_ant_MLGA= 378; %[K]
T_mean_ant_MLGA= 258; %[K]

% on this way, we are considering that T_equivalent of
% atmosphere/environment is : T_eq=(1-A)*T;
% with (A=0, so Max Attenuation in linear scale)

% Equivalent Noise Temperature of Transmission Lines (TL) and Waveguide(WG):
% ASSUME= A good TL and WG have a NF_TL= 1dB; (arbitrary)
NF_TL= 1; %[dB]
T0= 290; %[K]
T_TL=10^((NF_TL/10)-1)*T0; %[K]

% Equivalent Noise Temperature for LNA (Low Noise Amplifier):
% Typical NF for a good LNA is: NF_LNA=1dB
NF_LNA= 1; %[dB]
T0= 290; %[K]
T_LNA=10^((NF_LNA/10)-1)*T0; %[K]


%Total Temperature of the System :
%Case 1: Lander (Tx) - Mothership (Rx)
T_sys1a= T_mean_ant_HGA + T_TL + T_LNA ; %[K] Equivalent Noise Temperature of Overall System
T_sys1b= T_mean_ant_MLGA + T_TL + T_LNA; %[K] Equivalent Noise Temperature of Overall System
% T_sys1a= T_min_ant_HGA + T_TL + T_LNA ; %[K] BEST CASE Equivalent Noise Temperature of Overall System
% T_sys1b= T_min_ant_MLGA + T_TL + T_LNA; %[K] BEST CASE Equivalent Noise Temperature of Overall System

%Considering also the Noise due to the environment we add the temperature
%of surface of callisto, because the Rx is whatchin Callisto's surface
%without having the Sun in the Field of View :
T_eqN1a=T_sys1a+ T_env_rx1; %[K] Case 1 Uplink-Lander for HGA
T_eqN1b=T_sys1b+ T_env_rx1; %[K] Case 1 Uplink-Lander for MGA and LGA

%in case of DOWNLINK (Case 2) we have:

%Case 2: Mothership (Tx) - Lander (Rx)
%ASSUME: that the LNA of Lander is worst than the one on Mothership and we
%assume a CF=Correction Factor of CF=1.3 (arbitrary)
CF=1.3;
%operating values of lander antenna:
tmin1=233.15; %[K]
tmax1=333.15; %[K]
tmedio1=283.15; %[K]

T_sys2= tmedio1 + T_TL + T_LNA*CF ; %[K] Equivalent Noise Temperature of Overall System

%Considering also the environment noise, the Rx is whatching to the deep
%space, so we ASSUME that Sun and Jupiter are not present in the Field of
%View, or in the pointing direction/or in theta-3dB:
T_eqN2=T_sys2 + T_env_rx2; %[K] Case 2 Downlink-Lander 


%FINALLY we compute the NOISE for this 2 cases:
% B=0.01; %[GHz] (GHz=1e9Hz) Bandwidth of communication of X-band OMNIDIRECTIONAL for Uplink and Downlink
B=0.15; %[GHz] (GHz=1e9Hz) Bandwidth of communication of X-band for Uplink and Downlink
%Bandwidth of frequency range considered: 50 MHz(+-25 MHz from a center value), professor Luini sayd B=1GHz fro x-band.
Kb= 1.38e-23; %[W/Hz*K] Boltzmann Constant
Kb_dB= -228.6; %[dB/Hz*K] Boltzmann Constant 

%Noise power and dB, approximated values:
P_noise1a= Kb*T_eqN1a*(B*1e9); %[W] Case 1 Uplink-Lander for HGA
P_noise1a_dB= 10*log10(P_noise1a); %[dB]

P_noise1b= Kb*T_eqN1b*(B*1e9); %[W] Case 1 Uplink-Lander for MGA and LGA
P_noise1b_dB= 10*log10(P_noise1b); %[dB]

P_noise2= Kb*T_eqN2*(B*1e9); %[W]
P_noise2_dB= 10*log10(P_noise2); %[dB] Case 2 Downlink-Lander

%[From average case and Best case there are almost 5dB of Noise of difference !!!]

%Result:
figure()

subplot(1,3,1)
plot(r,P_rx1_dB(:,:),r,linspace(P_noise1b_dB,P_noise1b_dB,length(r)),'LineWidth',1.2)
title('Power_{Rx} [dB] for LGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('Recived Power [dB]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','P_{noise}= [dB]','Location','northeast')
grid on

subplot(1,3,2)
plot(r,P_rx2_dB(:,:),r,linspace(P_noise1b_dB,P_noise1b_dB,length(r)),'LineWidth',1.2)
title('Power_{Rx} [dB] for MGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('Recived Power [dB]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','P_{noise}= [dB]','Location','northeast')
grid on

subplot(1,3,3)
plot(r,P_rx3_dB(:,:),r,linspace(P_noise1a_dB,P_noise1a_dB,length(r)),'LineWidth',1.2)
title('Power_{Rx} [dB] for HGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('Recived Power [dB]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','P_{noise}= [dB]','Location','northeast')
grid on


figure ()
subplot(1,3,1)
plot(r,P_rx_1(:,:),r,linspace(P_noise1b,P_noise1b,length(r)),'LineWidth',1.2)
title('Power_{Rx} [W] for LGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('Recived Power [W]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','P_{noise}= [W]','Location','northeast')
grid on

subplot(1,3,2)
plot(r,P_rx_2(:,:),r,linspace(P_noise1b,P_noise1b,length(r)),'LineWidth',1.2)
title('Power_{Rx} [W] for MGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('Recived Power [W]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','P_{noise}= [W]','Location','northeast')
grid on

subplot(1,3,3)
plot(r,P_rx_3(:,:),r,linspace(P_noise1a,P_noise1a,length(r)),'LineWidth',1.2)
title('Power_{Rx} [W] for HGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('Recived Power [W]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','P_{noise}= [W]','Location','northeast')
grid on


%% LINK BUDGET COMPUTATION :

% Time of flyby : (Data from Europa Clipper Trajectory (pdf) - Trajectory =13F7-A21)
% Computed (by hand) for different flybys of Europa clipper: 1)h=25km; 2)h=466km; 3)h=1828km;
time_f1 = 68.8; %[s] h=25km
time_f2 = 730.8; %[s] h=466km
time_f3 = 2474.3; %[s] h=1828km
time_f= (time_f1+time_f2+time_f3)/3; %[s] average time of flybys

%Data quantity : (arbitrary)
Data_amount= 8e9; %[bit] equeal to 1Gbyte.
%Computation of Data Rate needed to deliver Data_amount in time_f of flyby coperture:
Data_rate= Data_amount/time_f1; %[bps]=bit per second
R=Data_rate*1.25; %[bit/s] + a safety factor of x1.25
% Data_rate is equal to =3.7 [Mbit/s] and the limit for our device is quite higher(150Mbps) - (from endurosat website)

%Computation Link Budget:
% P_rx3_dB(j,:)=P_tx_dB(j)+G_tx_dB+f_tx_dB+L_tx_dB+L_ptx_dB + (-L_s+L_mis+A_dB) + G_rx3_dB+f_rx_dB+L_rx_dB+L_prx3_dB;
% P_rx_1=10.^(P_rx1_dB(:,:)./10); %[W] Receiver_1 =LGA
% P_rx_2=10.^(P_rx2_dB(:,:)./10); %[W] Receiver_2 =MGA
% P_rx_3=10.^(P_rx3_dB(:,:)./10); %[W] Receiver_3 =HGA
% B=0.05; %[GHz] (GHz=1e9Hz) Bandwidth of communication of X-band for Uplink and Downlink
% Kb= 1.38e-23; %[W/Hz*K] Boltzmann Constant
% Kb_dB= -228.6; %[dB/Hz*K] Boltzmann Constant 
%Case 1: Lander (Tx) - Mothership (Rx)
% T_sys1a= T_mean_ant_HGA + T_TL + T_LNA ; %[K] Equivalent Noise Temperature of Overall System
% T_sys1b= T_mean_ant_MLGA + T_TL + T_LNA; %[K] Equivalent Noise Temperature of Overall System
% T_eqN1a=T_sys1a+ T_env_rx1; %[K] Case 1 Uplink-Lander for HGA
% T_eqN1b=T_sys1b+ T_env_rx1; %[K] Case 1 Uplink-Lander for MGA and LGA

%Computation of a normalized SNR (Signal To Noise Ratio [per bit]) for our antenna=
% used to choose properly the Modulation and Encoding alghoritms...
% (related to a particular graph, in books, and BER (Bit Error Ratio))

SNR1_bit= P_rx_1./(Kb*T_eqN1b*R); % [W] Receiver_1 =LGA
SNR1_dB_bit=10*log10(SNR1_bit);%[dB]
SNR2_bit= P_rx_2./(Kb*T_eqN1b*R); % [W] Receiver_2 =MGA
SNR2_dB_bit=10*log10(SNR2_bit);%[dB]
SNR3_bit= P_rx_3./(Kb*T_eqN1a*R); % [W] Receiver_3 =HGA
SNR3_dB_bit=10*log10(SNR3_bit);%[dB]
%THESE ARE THE VALUE OF OUR DISEGNED ANTENNA FOR A DATA RATE (R) OF 142 Mbit/s



%These values of SNR_dB_bit (or Normalized SNR) allow us to chose easily an
%Encoding and Modulation Algorithms ... (see slide to know it)

% Results
figure ()
subplot(1,3,1)
plot(r,SNR1_dB_bit(:,:),'LineWidth',1.2)
title('Normalized SNR (E_b/N_0) of trasmission - LGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('SNR per bit (E_b/N_0) [dB]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','Location','northeast')
grid on

subplot(1,3,2)
plot(r,SNR2_dB_bit(:,:),'LineWidth',1.2)
title('Normalized SNR (E_b/N_0) of trasmission - MGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('SNR per bit (E_b/N_0) [dB]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','Location','northeast')
grid on

subplot(1,3,3)
plot(r,SNR3_dB_bit(:,:),'LineWidth',1.2)
title('Normalized SNR (E_b/N_0) of trasmission - HGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('SNR per bit (E_b/N_0) [dB]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','Location','northeast')
grid on


%% Computation of effective SNR (Signal to Noise Ratio) :
%a bit useless...
SNR1= P_rx_1./(Kb*T_eqN1b*(B*1e9)); % [W] Receiver_1 =LGA
SNR1_dB=10*log10(SNR1);%[dB]
SNR2= P_rx_2./(Kb*T_eqN1b*(B*1e9)); % [W] Receiver_2 =MGA
SNR2_dB=10*log10(SNR2);%[dB]
SNR3= P_rx_3./(Kb*T_eqN1a*(B*1e9)); % [W] Receiver_3 =HGA
SNR3_dB=10*log10(SNR3);%[dB]

figure ()
subplot(1,3,1)
plot(r,SNR1_dB(:,:),'LineWidth',1.2)
title('SNR of trasmission - LGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('SNR [dB]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','Location','northeast')
grid on

subplot(1,3,2)
plot(r,SNR2_dB(:,:),'LineWidth',1.2)
title('SNR of trasmission - MGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('SNR [dB]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','Location','northeast')
grid on

subplot(1,3,3)
plot(r,SNR3_dB(:,:),'LineWidth',1.2)
title('SNR of trasmission - HGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('SNR [dB]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','Location','northeast')
grid on
 

%% Computation of Avaiable (maximum) Data rate for our modelled P_Rx:
%ASSUME A GOOD SNR: (data from internet wifi)
% no connectio :    SNR < 10 dB
% very low signal:  SNR = 10/15 dB
% low signal:       SNR = 15/25 dB
% very good signal: SNR = 25/40 dB
% excellent signal: SNR > 40 dB

SNR_fix_pb_dB=25; %[dB] - arbitrary
SNR_fix_per_bit= 10^(SNR_fix_pb_dB/10); %[W]

R_max1= P_rx_1./(Kb*T_eqN1b*SNR_fix_per_bit); %[bps]=bit per second Receiver_1 =LGA
R_max1_plot=R_max1*1e-6; %[Mbps]
R_max2= P_rx_2./(Kb*T_eqN1b*SNR_fix_per_bit); %[bps]=bit per second Receiver_2 =MGA
R_max2_plot=R_max2*1e-6; %[Gbps]
R_max3= P_rx_3./(Kb*T_eqN1a*SNR_fix_per_bit); %[bps]=bit per second Receiver_3 =HGA
R_max3_plot=R_max3*1e-9; %[Gbps]
%NOTE: our device limit is around 150Mbps.

%Results:
figure ()
subplot(1,3,1)
plot(r,R_max1_plot(:,:),'LineWidth',1.2)
title('Maximum Data rate for SNR=25dB, LGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('Data Rate [Mbps]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','Location','northeast')
grid on

subplot(1,3,2)
plot(r,R_max2_plot(:,:),'LineWidth',1.2)
title('Maximum Data rate for SNR=25dB, MGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('Data Rate [Mbps]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','Location','northeast')
grid on

subplot(1,3,3)
plot(r,R_max3_plot(:,:),'LineWidth',1.2)
title('Maximum Data rate for SNR=25dB, HGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('Data Rate [Gbps]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','Location','northeast')
grid on


%% HOW TO CHOSE A PROPER BANDWIDTH FOR COMMUNICATION

% For converting the ratio of two power values in dB, we use ans(dB) = 10*log10(ratio) and ratio=10^(ans(dB)/10).
% For converting the ratio of two voltage or current values in dB, we use ans(dB)=20*log10(ratio) and ratio=10^(ans(dB)/20). 
% 
% CR= 100e6; %[bit/s] (equal to 4Mbit/s, minimu value of Data transfer
% SN_dB= -20; %[dB] =SNR ipotetico
% SN=10^(SN_dB/10); %[W] =SNR ipotetico
% BB=CR/(log2(1+SN))   %[Hz] Bandwidth needed for this capacity channel (Shannon-Hartley Theorem)
% %resulting BB for SN=25dB and CR=4Mbit/s is: 5e5 [Hz], so: 0.0005 [GHz]=1e9
% %resulting BB for SN=25dB and CR=100Mbit/s is: 1.2e7 [Hz], so: 0.012 [GHz]=1e9

%Maximum Data Rate (from Shannon-Hartley Theorem):
MaxRate1=(B*1e9)*log2(1+SNR1); %[bits/s] =Channel Capacity
MaxRate2=(B*1e9)*log2(1+SNR2); %[bits/s] =Channel Capacity
MaxRate3=(B*1e9)*log2(1+SNR3); %[bits/s] =Channel Capacity
%output:
MaxCC1=max(max(MaxRate1))/1e6; %[Mbits/s]
MaxCC2=max(max(MaxRate2))/1e6; %[Mbits/s]
MaxCC3=max(max(MaxRate3))/1e6; %[Mbits/s]
MinCC1=min(min(MaxRate1))/1e6; %[Mbits/s]
MinCC2=min(min(MaxRate2))/1e6; %[Mbits/s]
MinCC3=min(min(MaxRate3))/1e6; %[Mbits/s]

%With the chosed value
Min_value_SNR_per_bit=2^((R/(B*1e9))-1)/(R/(B*1e9)); %[dB]
Channel_Capacity_our_antenna=[MaxCC1;MinCC1;MaxCC2;MinCC2;MaxCC3;MinCC3];
T=table(Channel_Capacity_our_antenna,'RowNames',{'Maximum CC [Mbit/s]- LGA','Minimum CC [Mbit/s]- LGA','Maximum CC [Mbit/s]- MGA','Minimum CC [Mbit/s]- MGA','Maximum CC [Mbit/s]- HGA','Minimum CC [Mbit/s]- HGA'})


%Results:
figure ()
subplot(1,3,1)
plot(r,MaxRate1(:,:)/1e3,'LineWidth',1.2)
title('Maximum Data rate, LGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('Maximum Channel Capacity [Kbps]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','Location','northeast')
grid on

subplot(1,3,2)
plot(r,MaxRate2(:,:)/1e6,'LineWidth',1.2)
title('Maximum Data rate, MGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('Maximum Channel Capacity [Mbps]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','Location','northeast')
grid on

subplot(1,3,3)
plot(r,MaxRate3(:,:)/1e6,'LineWidth',1.2)
title('Maximum Data Rate, HGA Receiver');
% title('Maximum Data Rate Omni-Trasmitter, HGA Receiver');
xlabel('Maximum travel of signal [km]');
ylabel('Maximum Channel Capacity [Mbps]');
legend('P_{tx}=1W','P_{tx}=2W','P_{tx}=3W','P_{tx}=4W','Location','northeast')
grid on



%% TRIAL, MAL RIUSCITI...

% SNR1_bit= P_rx_1./(Kb*T_eqN1b*R); % [W] Receiver_1 =LGA
% SNR1_dB_bit=10*log10(SNR1_bit);%[dB]
% SNR2_bit= P_rx_2./(Kb*T_eqN1b*R); % [W] Receiver_2 =MGA
% SNR2_dB_bit=10*log10(SNR2_bit);%[dB]
% SNR3_bit= P_rx_3./(Kb*T_eqN1a*R); % [W] Receiver_3 =HGA
% SNR3_dB_bit=10*log10(SNR3_bit);%[dB]
%THESE ARE THE VALUE OF OUR DISEGNED ANTENNA FOR A DATA RATE (R) OF 3.7Mbit/s
% 
% fix_R= fix_B; %[bit/s] = 5 Mbit/s

% myfun= @(x,u) (log2(x)./x)-u;
% BB1=(1+SNR11)/fix_R;
% BB2=(1+SNR22)/fix_R;
% BB3=(1+SNR33)/fix_R;

% [low_B1,FVAL,EXITFLAG,OUTPUT]=fsolve(@(x) myfun(x,BB1),10000);
% [low_B2,FVAL,EXITFLAG,OUTPUT]=fsolve(@(x) myfun(x,BB2),10000);
% [low_B3,FVAL,EXITFLAG,OUTPUT]=fsolve(@(x) myfun(x,BB3),100000);
% 
% x=[1e6:1000:1e9];
% y=(log2(x)./x)-BB1;
% figure
% plot(x,y)

