function [x1, y1, x2, y2] = colouredplot(lon, lat, R, h0, FOV, i0)

% hold all
% plot(x,y1)
% plot(x,y2)
% patch([x fliplr(x)], [y1 fliplr(y2)], 'g')
% hold off


rho = asin(R/(R+h0));
% lambda_0 = acos(R/(R+h0));

eps = acos(sin(FOV)/sin(rho));
lambda = pi/2-FOV-eps;

x1 = zeros(length(lon),1);
x2 = zeros(length(lon),1);
y1 = zeros(length(lon),1);
y2 = zeros(length(lon),1);


for ii= 1:length(lon)
x1(ii) = -sign(sin(lat(ii)+(pi-i0))).*(rad2deg(lambda))*abs(cos(pi/2*deg2rad(lat(ii))/i0))+ lon(ii);
x2(ii) = +sign(sin(lat(ii)+(pi-i0))).*(rad2deg(lambda))*abs(cos(pi/2*deg2rad(lat(ii))/i0)) + lon(ii);
y1(ii) = +rad2deg(lambda)*abs(sin(pi/2*deg2rad(lat(ii))/i0)) + lat(ii);
y2(ii) = -rad2deg(lambda)*abs(sin(pi/2*deg2rad(lat(ii))/i0)) + lat(ii);
ii = ii +1;
end



figure
hold all
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon,lat,'g','linewidth',1.2);
plot(lon(1),lat(1),'go','linewidth',2)
plot(lon(length(lon)),lat(length(lat)),'gs','linewidth',2)
% s = patch([x1 fliplr(x1)], [y1 fliplr(y2)], 'b');
plot(x1,y1,'b','linewidth',1.2)
plot(x2,y2,'b','linewidth',1.2)
xlim([-180 180])
ylim([-90 90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('Ground track','Start','End', 'Field of View', 'Location','best','Orientation','horizontal')
title('Ground track plot')
s.FaceVertexAlphaData = 1;
s.FaceAlpha = 'flat' ; 
hold off

end