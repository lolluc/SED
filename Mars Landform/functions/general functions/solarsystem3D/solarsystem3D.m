% multi-purpose 3D planet plot generator for the entire solar system
%
%
% -----------INPUTS----------
% planet_num= 0 SUN
%             1 MERCURY
%             2 VENUS
%             3 EARTH  (default if left blank)
%             3.5 MOON
%             4 MARS
%             5 JUPITER
%             6 SATURN
%             7 URANUS
%             8 NEPTUNE
%
% color= 'black' - uses milky way background and white axis and labels
%        'white' - no backgroung
%
% r_sun= vector of the sun seen from the planet (used for lighting)
%
% -----------OUTPUTS---------
% 3D planet plot with labelled axis
%
%

function []= solarsystem3D(planet_num,color,r_sun)

%settings------ (automatically uses earth if there are no inputs)
if nargin== 0
    planet_num=3;
    color='black';
    r_sun=[0 0 0];
end
if nargin== 1
    color='black';
    r_sun=[0 0 0];
end
if nargin== 2
    r_sun=[0 0 0];
end
npanels = 300;% number of panels of planets spheres

%figure making--------------------------
figure
hold on;
grid on;

% Set the axes scale equal
axis equal;

% Put the axes labels
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
% Set initial view angle
view(120,30);

%background figure (milky way) setting--------------
switch color
    case 'black'
        % This creates the 'background' axes
        ha = axes('units','normalized','position',[0 0 1 1]);
        % Move the background axes to the bottom
        uistack(ha,'bottom');
        % Load in a background image and display it using the correct colors
        I=imread('sky.jpg'); imagesc(I);
        % Turn the handlevisibility off so that we don't inadvertently plot into
        % the axes again. Also, make the axes invisible
        set(ha,'handlevisibility','off', 'visible','off')
end
% various planet setting and texturing----------------------


% PLANET 0: SUN  +++++++++++++
if planet_num==0
    R=696340*10; %sun radius
    [x,y,z]= ellipsoid(0, 0, 0, R, R, R, npanels);
    planet = surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0.04,'DiffuseStrength',1,'HandleVisibility','off');
    % Load sun image for texture map
    cdata = imread('sun.jpg');
    % Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
    % specify the image data using the 'CData' property with the data loaded
    % from the image. Finally, set the transparency and remove the edges.
    set(planet, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    shading interp
    
    % PLANET 1: MERCURY  +++++++++++++++++
    
    
elseif planet_num==1
    R=2439.7; %mercury radius
    [x,y,z]= ellipsoid(0, 0, 0, R, R, R, npanels);
    planet = surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0.04,'DiffuseStrength',1,'HandleVisibility','off');
    % Load mercury image for texture map
    cdata = imread('mercury.jpg');
    % Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
    % specify the image data using the 'CData' property with the data loaded
    % from the image. Finally, set the transparency and remove the edges.
    set(planet, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    shading interp
    
    %lightning
    if norm(r_sun)==0
        lux=[0,cosd(1/30)*62.869*10^6,sind(1/30)*62.869*10^6]; %light position
    else
        lux=r_sun; %light position
    end
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    
    
    
    
    % PLANET 2: VENUS  +++++++++++++++++
    
elseif planet_num==2
    R=6051.8; %venus radius
    [x,y,z]= ellipsoid(0, 0, 0, R, R, R, npanels);
    planet = surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0.04,'DiffuseStrength',1,'HandleVisibility','off');
    % Load venus image for texture map
    cdata = imread('venus_surf.jpg');
    % Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
    % specify the image data using the 'CData' property with the data loaded
    % from the image. Finally, set the transparency and remove the edges.
    set(planet, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    shading interp
    
    % cloud texture
    cdata2=imread('venus_atm.jpg');
    
    %cloud surface sphere (100 km heigth)
    [xatm, yatm, zatm] = ellipsoid(0, 0, 0, R+100, R+100, R+100, npanels);
    atmo = surf(xatm, yatm, -zatm, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0.01,'DiffuseStrength',1,'HandleVisibility','off');
    %apply clouds texture and set transparency
    set(atmo, 'FaceColor', 'texturemap', 'CData', cdata2 ,...
        'FaceAlpha',0.85,'EdgeColor', 'none');
    shading interp
    
    %lightning
    if norm(r_sun)==0
        lux=[0,cosd(3)*107.63*10^6,sind(3)*107.63*10^6]; %light position
    else
        lux=r_sun; %light position
    end
    
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    
    
    
    
    % PLANET 3: EARTH  +++++++++++++++++
    
elseif planet_num==3
    R=6371; %earth radius
    [x,y,z]= ellipsoid(0, 0, 0, R, R, R, npanels);
    planet = surf(x,y,-z, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0.04,'DiffuseStrength',1,'HandleVisibility','off');
    % Load Earth image for texture map
    cdata = imread('earth.tif');
    % Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
    % specify the image data using the 'CData' property with the data loaded
    % from the image. Finally, set the transparency and remove the edges.
    set(planet, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    shading interp
    
    
    % clouds texture and transparency
    alpha2=importdata('clouds.png');
    cdata2=imread('clouds.png');
    alpha2=alpha2.alpha;
    alpha2=double(alpha2)/255;
    %clouds surface sphere (100 km heigth)
    [xatm, yatm, zatm] = ellipsoid(0, 0, 0, R+100, R+100, R+100, npanels);
    atmo = surf(xatm, yatm, -zatm, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0,'DiffuseStrength',1,'HandleVisibility','off');
    %apply clouds texture
    set(atmo, 'FaceColor', 'texturemap', 'CData', cdata2 ,...
        'FaceAlpha','texturemap','AlphaData',alpha2);
    shading interp
    
    %blue atmosphere filter sphere (300 km heigth)
    [xfil, yfil, zfil] = ellipsoid(0, 0, 0, R+300, R+300, R+300, npanels);
    filter = surf(xfil, yfil, -zfil, 'FaceColor', [37 60 190]/255, 'EdgeColor',...
        'none','SpecularStrength',0,'DiffuseStrength',0.5,'HandleVisibility','off');
    set(filter,'FaceAlpha',0.2,'AmbientStrength',0.4);
    
    
    %lightning
    if norm(r_sun)==0
        lux=[0,cosd(23.5)*150*10^6,sind(23.5)*150*10^6]; %light position
    else
        lux=r_sun; %light position
    end
    
    light('position',lux,'style','local',...
        'color',[255 255 233]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 233]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 233]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 233]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 233]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 233]/255,'HandleVisibility','off');
    
    
    
    
    % PLANET 3.5: MOON  +++++++++++++++++
    
elseif planet_num==3.5
    R=1737.1; %moon radius
    [x,y,z]= ellipsoid(0, 0, 0, R, R, R, npanels);
    planet = surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0.04,'DiffuseStrength',1,'HandleVisibility','off');
    % Load moon image for texture map
    cdata = imread('moon.jpg');
    % Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
    % specify the image data using the 'CData' property with the data loaded
    % from the image. Finally, set the transparency and remove the edges.
    set(planet, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    shading interp
    
    %lightning
    if norm(r_sun)==0
        lux=[0,cosd(28)*150*10^6,sind(28)*150*10^6]; %light position
    else
        lux=r_sun; %light position
    end
    
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    
    
    
    % PLANET 4: MARS  +++++++++++++++++
    
elseif planet_num==4
    R=3389.5; %mars radius
    [x,y,z]= ellipsoid(0, 0, 0, R, R, R, npanels);
    planet = surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0.02,'DiffuseStrength',1,'HandleVisibility','off');
    % Load mars image for texture map
    cdata = imread('MarsTexture.jpg');
    % Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
    % specify the image data using the 'CData' property with the data loaded
    % from the image. Finally, set the transparency and remove the edges.
    set(planet, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    shading interp
    
    %orange-ish atmosphere filter sphere (50 km heigth)
    [xfil, yfil, zfil] = ellipsoid(0, 0, 0, R+50, R+50, R+50, npanels);
    filter = surf(xfil, yfil, -zfil, 'FaceColor', [210 226 251]/255, 'EdgeColor',...
        'none','SpecularStrength',0,'DiffuseStrength',0.5,'HandleVisibility','off');
    set(filter,'FaceAlpha',0.2,'AmbientStrength',0.4);
    
    %lightning
    if norm(r_sun)==0
        lux=[0,cosd(25)*212.3*10^6,sind(25)*212.3*10^6]; %light position
    else
        lux=r_sun; %light position
    end
    
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    
    
    
    
    % PLANET 5: JUPITER  +++++++++++++++++
    
elseif planet_num==5
    R=69911; %jupiter radius
    [x,y,z]= ellipsoid(0, 0, 0, R, R, R, npanels);
    planet = surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0.02,'DiffuseStrength',1,'HandleVisibility','off');
    % Load jupiter image for texture map
    cdata = imread('jupiter.jpg');
    % Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
    % specify the image data using the 'CData' property with the data loaded
    % from the image. Finally, set the transparency and remove the edges.
    set(planet, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    shading interp
    
    %lightning
    if norm(r_sun)==0
        lux=[0,cosd(3)*767*10^6,sind(3)*767*10^6]; %light position
    else
        lux=r_sun; %light position
    end
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    
    
    
    
    % PLANET 5: SATURN  +++++++++++++++++
    
elseif planet_num==6
    R=58232; %saturn radius
    [x,y,z]= ellipsoid(0, 0, 0, R, R, R, npanels);
    planet = surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0.02,'DiffuseStrength',1,'HandleVisibility','off');
    % Load saturn image for texture map
    cdata = imread('saturn_surf.jpg');
    % Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
    % specify the image data using the 'CData' property with the data loaded
    % from the image. Finally, set the transparency and remove the edges.
    set(planet, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    shading interp
    
    % ring texture and transparency
    alpha2=importdata('saturn_ring.png');
    cdata2=imread('saturn_ring.png');
    alpha2=alpha2.alpha;
    alpha2=double(alpha2)/255;
    angle=0:0.02*pi:2*pi;
    xring = zeros(2,length(angle)); yring= xring;
    for i=1:length(angle)
        xring(:,i) =  ([R,140220] * cos(angle(i)))' ;
        yring(:,i) =  ([R,140220]*  sin(angle(i)))' ;
        
    end
    zring = zeros(size(xring));
    ring = surf(xring, yring, zring, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0.1,'DiffuseStrength',1,'HandleVisibility','off');
    %apply rings texture
    set(ring, 'FaceColor', 'texturemap', 'CData', cdata2 ,...
        'FaceAlpha','texturemap','AlphaData',alpha2,'ambientstrength',1);
    shading flat
    
    
    
    %lightning
    if norm(r_sun)==0
        lux=[0,cosd(26.3)*1493*10^6,sind(26.3)*1493*10^6]; %light position
    else
        lux=r_sun; %light position
    end
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    
    
    
    
    % PLANET 7: URANUS  +++++++++++++++++
    
elseif planet_num==7
    R=25362; %uranus radius
    [x,y,z]= ellipsoid(0, 0, 0, R, R, R, npanels);
    planet = surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0.04,'DiffuseStrength',1,'HandleVisibility','off');
    % Load uranus image for texture map
    cdata = imread('uranus.jpg');
    % Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
    % specify the image data using the 'CData' property with the data loaded
    % from the image. Finally, set the transparency and remove the edges.
    set(planet, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    shading interp
    
    %lightning
    if norm(r_sun)==0
        lux=[0,cosd(98)*2959*10^6,sind(98)*2959*10^6]; %light position
    else
        lux=r_sun; %light position
    end
    
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    
    
    
    
    % PLANET 8: NEPTUNE  +++++++++++++++++
    
elseif planet_num==8
    R=24622; %neptune radius
    [x,y,z]= ellipsoid(0, 0, 0, R, R, R, npanels);
    planet = surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 'none',...
        'SpecularStrength',0.04,'DiffuseStrength',1,'HandleVisibility','off');
    % Load neptune image for texture map
    cdata = imread('neptune.jpg');
    % Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
    % specify the image data using the 'CData' property with the data loaded
    % from the image. Finally, set the transparency and remove the edges.
    set(planet, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    shading interp
    
    %lightning
    if norm(r_sun)==0
        lux=[0,cosd(29)*4476*10^6,sind(29)*4476*10^6]; %light position
    else
        lux=r_sun; %light position
    end
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    light('position',lux,'style','local',...
        'color',[255 255 253]/255,'HandleVisibility','off');
    
    
    
else   %in case the input planet number is set wrong:
    error('wrong planet number!')
end


% axes color and transparency-------------
switch color
    case 'black'
        ax = gca;
        ax.XColor = [200 200 200]/255;
        ax.YColor = [200 200 200]/255;
        ax.ZColor = [200 200 200]/255;
        ax.GridColor = 'white';
        ax.GridAlpha = 0.27;
        ax.Color='none';
end
return

