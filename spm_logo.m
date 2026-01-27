function h = spm_logo(S)
% Generate SPM's new logo
% FORMAT h = spm_logo(S)
%
% S           - input structure
%  Optional fields of S:
%   S.ver     - Set the version number, defaults to current SPM version
%   S.altm    - Which style of M to use? False for M and True for m.
%   S.width   - 2x1 array to set widths of the text of logo.
%   S.colour  - Text colour
%   S.background    - background colour
%   S.angle   - Angle (in degrees) of the isometric progection
%   S.save    - If true, generates image files for re-use
%
% h           - Handle figure
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2025 UCL Dept. of Imaging Neuroscience

if ~nargin
    S = [];
end

if ~isfield(S,'ver'), S.ver = spm('ver'); end
if ~isfield(S,'altm'), S.altm = true; end
if ~isfield(S,'width'), S.width = [0.05 0.02]; end
if ~isfield(S,'colour'), S.colour = 'k'; end
if ~isfield(S,'background'), S.background = 'w'; end
if ~isfield(S,'angle'), S.angle = -30; end
if ~isfield(S,'save'),  S.save = false; end

if isscalar(S.width), S.width = [S.width S.width]; end
if isnumeric(S.ver), S.ver = num2str(S.ver); end
if strcmp(S.ver(1),'S'), S.ver = S.ver(4:end); end

% Define digits (0–9) as line segments
digits = {
    'S', [0 0; 1 0; 1 0.5; 0 0.5; 0 1; 1 1];
    'P', [0 0; 0 1; 0 1; 1 1; 1 0.5; 0 0.5];
    'M', [0 0; 0 1; 0.5 0.5; 1 1; 1 0];
    'm', [0 0; 0 1; 0.5 1; 0.5 0; 0.5 1; 1 1; 1 0];
    '0', [0 0; 0 1; 1 1; 1 0; 0 0];
    '1', [1 0; 1 1];
    '2', [0 1; 1 1; 1 0.5; 0 0.5; 0 0; 1 0];
    '3', [0 1; 1 1; 1 0.5; 0 0.5; 1 0.5; 1 0; 0 0];
    '4', [1 0; 1 1; 1 0.5; 0 0.5; 0 1];
    '5', [1 1; 0 1; 0 0.5; 1 0.5; 1 0; 0 0];
    '6', [1 1; 0 1; 0 0; 1 0; 1 0.5; 0 0.5];
    '7', [0 1; 1 1; 0.5 0];
    '8', [0 0; 0 1; 1 1; 1 0; 0 0; 0 0.5; 1 0.5];
    '9', [0 0; 1 0; 1 1; 0 1; 0 0.5; 1 0.5;];
    };


h = figure;
hold on;
axis equal off;
offset = 0;

if S.altm
    digit_string = 'SPm';
else
    digit_string = 'SPM';
end

for l = 1:length(digit_string)
    digit = digit_string(l);
    coords = digits{strcmp(digit, digits(:,1)), 2};  % Find the digit coordinates

    % Plot each pair of coordinates as a line segment
    for i = 1:size(coords, 1)-1
        x0 = coords(i, 1) + offset;
        z0 = coords(i, 2);
        x1 = coords(i + 1, 1) + offset;
        z1 = coords(i + 1, 2);
        y = 0;  % constant depth
        [ix0, iz0] = isometric_transform_3d(x0, y, z0, S.angle);
        [ix1, iz1] = isometric_transform_3d(x1, y, z1, S.angle);
        c = line2rec([ix0 iz0],[ix1 iz1],S.width(1));
        fill(c(:,1),c(:,2),S.colour)
        plot_circle(ix0,iz0,S.width(1))
    end
    plot_circle(ix1,iz1,S.width(1))

    offset = offset + 1.25;  % Space between letters
end

set(h,'Color',S.background)

offset = 2.55;

digit_string = S.ver;

for l = 1:length(digit_string)
    digit = digit_string(l);
    coords = digits{strcmp(digit, digits(:,1)), 2};  % Find the digit coordinates

    % Plot each pair of coordinates as a line segment
    for i = 1:size(coords, 1)-1
        x0 = 0.4 * coords(i, 1) + offset;
        y0 = 0.6 * coords(i, 2) - 0.8;
        x1 = 0.4 * coords(i + 1, 1) + offset;
        y1 = 0.6 * coords(i + 1, 2) - 0.8;
        z = 0;  % constant depth
        [ix0, iy0] = isometric_transform_3d(x0, -y0, z, S.angle);
        [ix1, iy1] = isometric_transform_3d(x1, -y1, z, S.angle);
        c = line2rec([ix0 iy0],[ix1 iy1],S.width(2));
        fill(c(:,1),c(:,2),'k')
        plot_circle(ix0,iy0,S.width(2))
    end
    plot_circle(ix1,iy1,S.width(2))


    offset = offset + 0.6;  % Space between letters
end

% Write image(s) if requested
if S.save
    saveas(h,'spm_logo.png');
end

end

% isometric projection function
function [iso_x, iso_y] = isometric_transform_3d(x, y, z, ang)
if nargin < 4
    ang = -30;
end
angle = deg2rad(ang);
iso_x = x * cos(angle) - y * cos(angle);
iso_y = x * sin(angle) + y * sin(angle) + z;
end

% turn the line into a ractangle of some width
function corners = line2rec(p1,p2,w)
del = p2-p1;
% need width in orthogonal travel direction
orthdir =[-del(2), del(1)];
orthdir = orthdir./norm(orthdir);
corners = [p1-w*orthdir; p1+w*orthdir; p2+w*orthdir; p2-w*orthdir];
end

function plot_circle(x,y,r)
theta = linspace(0, 2*pi, 25);
x_circle = r * cos(theta);
y_circle = r * sin(theta);
fill(x + x_circle, y + y_circle, 'k', 'EdgeColor', 'none');
end