function [Rxo,Ryo,Rzo,Rxp,Ryp,Rzp] = plot_orbit(a,e,i,OM,om,vect,theta)
%Function che prende in input i parametri orbitali dell'orbita e restituisce i vettori delle sue componenti.
%Restituisce inoltre i valori per un punto dell'orbita.
%a = semiasse maggiore
%e = eccentricita'
%i = inclinazione in radianti
%OM = ascensione retta in radianti
%om = anomalia del pericentro in radianti
%theta = anomali punto sull'orbita in radianti

%% Scrivo equazione orbita
if vect == 0
alpha = linspace (0,2*pi,1000);  %Definisco un vettore alpha che descriva interamente l'orbita
else
alpha = vect;
end

p = (1-e^2) * a ;                                         %Definisco p per l'orbita 
erre = @(alpha) p ./ ( 1 + e*cos(alpha) );                %Scrivo l'equazione dell'orbita (ellisse) rispetto al sistema di riferimento perifocale
R= [erre(-alpha)];                                         %Definisco un vettore che descrive l'orbita calcolando l'equazione nei punti alpha

x= [R.*cos(alpha)];                                       %Scompongo R nelle componenti rispetto al sistema perifocale
y=[R.*sin(alpha)];
z= [R.*zeros];

R_p= [x; y; z];                                           %Definisco la matrice delle componenti dell'orbita nel sistem di riferimento perifocale

%% Trasformazione da sistema perifocale a sistema geocentrico orbita 
[T_gepe,T_pege] = rot(OM,om,i);

R_g= T_pege*R_p;


%% Plot orbita 
Rxo=R_g(1,:);
Ryo=R_g(2,:);
Rzo=R_g(3,:);


%% Plot punto 
if nargin == 7;

P= [erre(theta)];

xp= [P.*cos(theta)];
yp=[P.*sin(theta)];
zp= [P.*zeros];
PP = [xp;yp;zp];

P_g = T_pege*PP;

Rxp=P_g(1);
Ryp=P_g(2);
Rzp=P_g(3);

end


end