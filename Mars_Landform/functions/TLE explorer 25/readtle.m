function [a,e,i,OM,om,th,epoch2000]=readtle(file)

% READTLE Read satellite ephemeris data from a NORAD two-line element (TLE)
% file.
%
% INPUTS:
%   file    - Path to any standard two-line element file.
%
% OUTPUTS:
%   a[1xN]   - semi majour axis [km]
%   e[1xN]   - eccentricity [-]
%   i[1xN]   - inclination [deg]
%   OM[1xN]  - RAAN [deg]
%   om[1xN]  - argument of periapsis [deg]
%   th[1xN]  - real anomaly[deg]

mu=astroConstants(13);

fd = fopen(file,'r');
if fd < 0, fd = fopen([file '.tle'],'r'); end
assert(fd > 0,['Can''t open file ' file ' for reading.'])

n = 0;

A1 = fgetl(fd);
A2 = fgetl(fd);

while ischar(A2)
    n = n + 1;
    
    
    assert(chksum(A1), 'Checksum failure on line 1')
    assert(chksum(A2), 'Checksum failure on line 2')
    year= str2num(A1(19:20));
    if year<52
        year=year+2000;
    else
        year=year+1900;
    end
    epoch2000(n)=date2mjd2000([year,0,0,0,0,0])+str2num(A1(21:32));
    
    
    
    i(n) = str2num(A2(9:16));
    OM(n) = str2num(A2(18:25));
    e(n) = str2num(['.' A2(27:33)]);
    om(n) = str2num(A2(35:42));
    
    M = str2num(A2(44:51));
    E=fzero(@(E) -M+E-e(n)*sin(E),M);
    th(n)=2*atan(sqrt((1+e(n))/(1-e(n)))*tan(E/2));
    
    nn = str2num(A2(53:63));%rev/day
    T=86400/nn; %sec
    a(n) = ((T/(2*pi))^2*mu)^(1/3); %km
    
    
    A1 = fgetl(fd);
    A2 = fgetl(fd);
end

fclose(fd);

[epoch2000,II]=sort(epoch2000);
a=a(II);
e=e(II);
i=i(II);
OM=OM(II);
om=om(II);
th=th(II);


ii=1;
while ii<=length(a)-1
    ii=ii+1;
    if epoch2000(ii)==epoch2000(ii-1)
        epoch2000=[epoch2000(1:ii-1),epoch2000(ii+1:end)];
        a=[a(1:ii-1),a(ii+1:end)];
        e=[e(1:ii-1),e(ii+1:end)];
        i=[i(1:ii-1),i(ii+1:end)];
        OM=[OM(1:ii-1),OM(ii+1:end)];        
        om=[om(1:ii-1),om(ii+1:end)];
        th=[th(1:ii-1),th(ii+1:end)];
        ii=ii-1;
    end

end


th=rad2deg(th);
for ii=1:length(th)
    if th(ii)<0
        th(ii)=360+th(ii);
    end
end


end
%%
% Checksum (Modulo 10)
% Letters, blanks, periods, plus signs = 0; minus signs = 1
function result = chksum(str)
result = false; c = 0;

for k = 1:68
    if str(k) > '0' && str(k) <= '9'
        c = c + str(k) - 48;
    elseif str(k) == '-'
        c = c + 1;
    end
end

if mod(c,10) == str(69) - 48
    result = true;
end

end
