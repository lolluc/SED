function [Rxp,Ryp,Rzp] = plot_point(a,e,i,OM,om,theta)
%%Function che prende in input i parametri orbitali dell'orbita e 
%restituisce i valori del punto richiesto.
%a = semiasse maggiore
%e = eccentricita'
%i = inclinazione in radianti
%OM = ascensione retta in radianti
%om = anomalia del pericentro in radianti
%theta = anomali punto sull'orbita in radianti
%% Scrivo equazione orbita

alpha = theta;

p = (1-e^2) * a ;                                         %Definisco p per l'orbita 

P= p ./ ( 1 + e*cos(alpha) );
[T_gepe,T_pege] = rot(OM,om,i);

xp= [P.*cos(theta)];
yp=[P.*sin(theta)];
zp= [P.*zeros];
PP = [xp;yp;zp];

P_g = T_pege*PP;

Rxp=P_g(1);
Ryp=P_g(2);
Rzp=P_g(3);