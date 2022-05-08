function globe=Celestial_body(name,scale)

% Celestial_body.m - surface plot of a celestial body selected by the user.
%
% PROTOTYPE:
% globe=Celestial_body(name,scale).
%
% DESCRIPTION:
% Function to load the selected celestial body as a sphere inside a figure.
%
% INPUT:
% name  [string] name of the desired celestial body.
% scale [1x1]    scaling factor (optional: set to 1 by default)
% -available bodies:
%         sun
%         mercury
%         venus
%         earth
%         moon
%         mars
%         jupiter
%         saturn
%         uranus
%         neptune
%
% OUTPUT:
% globe [surface] surface of the selected celestial body
%
%   Notes for upgrading this function:
%       It is possible to change some operations to make it faster.
%       - DO NOT change the structure of the function, as well as its
%           prototype.
%       - DO NOT change the output type.
%       - DO NOT change the input to be given in degrees
%       Contact the authors for modification
%
% CALLED FUNCTIONS:
% (none)
%
% REFERENCE:
% radii taken from: Curtis, Howard D. Orbital Mechanics for Engineering Students. Butterworth-Heinemann, 2019.
%
% CHANGELOG:
% 23/11/2020, Riccardo Majer: Updated
% ------------------------------------------------------------------------------------------------------------

%add the textures folder to the current path

addpath(genpath('textures')); %NOTE: change the directory if textures folder and Celestial_body.m aren't in the same folder

% switch to select the body selected by the user
switch name
    case 'sun'
        R=696000;
        if nargin==2
            R = R*scale;
        end
        image =imread('sun.jpg');
    case 'mercury'
        R=2440;
        if nargin==2
            R=R*scale;
        end
        image =imread('mercury.jpg');
    case 'venus'
        R=6052;
        if nargin==2
            R=R*scale;
        end
        image =imread('venus.jpg');
    case 'earth'
        R=6378;
        if nargin==2
            R=R*scale;
        end
        image =imread('EarthTexture.jpg');
    case 'mars'
        R=3396;
        if nargin==2
            R=R*scale;
        end
        image =imread('MarsTexture.jpg');
    case 'jupiter'
        R=71490;
        if nargin==2
            R=R*scale;
        end
        image =imread('jupiter.jpg');
    case 'saturn'
        R=60270;
        if nargin==2
            R=R*scale;
        end
        image =imread('saturn.jpg');
    case 'uranus'
        R=25560;
        if nargin==2
            R=R*scale;
        end
        image =imread('uranus.jpg');
    case 'neptune'
        R=24764;
        if nargin==2
            R=R*scale;
        end
        image =imread('neptune.jpg');
    case 'moon'
        R=1737;
        if nargin==2
            R=R*scale;
        end
        image =imread('moon.jpg');
    otherwise
        error('%s is unavailable',name)
end

% Define the number of panels to be used to model the sphere
npanels=180;

% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x,y,z]=ellipsoid(0,0,0,R,R,R,npanels);

% Create the globe with the surf function
globe=surf(x,y,-z,'FaceColor','none','EdgeColor','none');

% Load celestial body image for texture map
cdata=(image);

% Set the transparency of the globe: 1 = opaque, 0 = invisible
alpha=1;

% Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
% specify the image data using the 'CData' property with the data loaded
% from the image. Finally, set the transparency and remove the edges
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none','AlignVertexCenters','on');