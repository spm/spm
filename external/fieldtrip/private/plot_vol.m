function plot_vol(vol, varargin)

% PLOT_VOL visualizes the boundaries in the vol structure constituting the
% geometrical information of the forward model
%
% Use as
%   hs = plot_vol(vol, varargin)
%
% Graphic facilities are available for vertices, edges and faces. A list of
% the arguments is given below with the correspondent admitted choices.
%
%     'faces'         true or false
%     'vertices'      true or false
%     'edges'         true or false
%     'facecolor'     [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'vertexcolor'   [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'edgecolor'     [r g b] values or string, for example 'brain', 'cortex', 'skin', 'black', 'red', 'r'
%     'faceindex'     true or false
%     'vertexindex'   true or false
%     'colormap'      'gray', 'hot', 'cool', 'copper', 'spring', 'summer', ...
%
% Example
%   vol.r = [86 88 92 100];
%   vol.o = [0 0 40];
%   figure, plot_vol(vol,'colormap','cool')

% Copyright (C) 2009, Cristiano Micheli
%
% $Log: plot_vol.m,v $
% Revision 1.8  2009/06/03 10:05:40  roboos
% changed handling of mesh generation and teh actual plotting
%
% Revision 1.7  2009/04/22 11:45:02  crimic
% updated help
%
% Revision 1.6  2009/04/22 11:43:57  crimic
% fixed colormap feature
%
% Revision 1.5  2009/04/19 11:35:16  crimic
% added colormap for all vol options
%
% Revision 1.4  2009/04/17 17:00:59  crimic
% updated help
%
% Revision 1.3  2009/04/17 16:48:34  crimic
% updated help
%
% Revision 1.2  2009/04/17 16:45:40  crimic
% created function to plot forward model geometry
%

% get the optional input arguments
faces       = keyval('faces',         varargin); if isempty(faces),       faces = true;       end
vertices    = keyval('vertices',      varargin); if isempty(vertices),    vertices = false;   end
edges       = keyval('edges',         varargin); if isempty(edges),       edges = true;      end
facecolor   = keyval('facecolor',     varargin); if isempty(facecolor),   facecolor = [0.5 0.5 0.5]; end % grey
facealpha   = keyval('facealpha',     varargin); if isempty(facealpha),   facealpha = 0.5;    end % slightly opaque
faceindex   = keyval('faceindex',     varargin);
vertexcolor = keyval('vertexcolor',   varargin);
vertexindex = keyval('vertexindex',   varargin);
vertexsize  = keyval('vertexsize',    varargin); if isempty(vertexsize),  vertexsize = 10;    end
edgecolor   = keyval('edgecolor',     varargin);
map         = keyval('colormap',      varargin);

% we will probably need a sphere, so let's prepare one
[pnt, tri] = icosahedron642;


% prepare a single or multiple triangulated boundaries
switch voltype(vol)
  case {'singlesphere' 'concentric'}
    vol.r = sort(vol.r);
    bnd = [];
    for i=1:length(vol.r)
      bnd(i).pnt(:,1) = pnt(:,1)*vol.r(i) + vol.o(1);
      bnd(i).pnt(:,2) = pnt(:,2)*vol.r(i) + vol.o(2);
      bnd(i).pnt(:,3) = pnt(:,3)*vol.r(i) + vol.o(3);
      bnd(i).tri = tri;
    end

  case 'multisphere'
    bnd = [];
    for i=1:length(vol.label)
      bnd(i).pnt(:,1) = pnt(:,1)*vol.r(i) + vol.o(i,1);
      bnd(i).pnt(:,2) = pnt(:,2)*vol.r(i) + vol.o(i,2);
      bnd(i).pnt(:,3) = pnt(:,3)*vol.r(i) + vol.o(i,3);
      bnd(i).tri = tri;
    end

  case {'bem', 'dipoli', 'asa', 'avo', 'bemcp', 'nolte'}
    % these already contain one or multiple triangulated surfaces for the boundaries
    bnd = vol.bnd;

  otherwise
    error('unsupported voltype')
end

if ~isempty(map)
  try
    cmap=colormap(map);close(gcf);
  catch
    error('Colormap does not exist')
  end
  % set the color of the sphere:
  r = (linspace(cmap(1,1),cmap(end,1),length(bnd)))';
  g = (linspace(cmap(1,2),cmap(end,2),length(bnd)))';
  b = (linspace(cmap(1,3),cmap(end,3),length(bnd)))';
  color = [r g b];
end

% plot the triangulated surfaces of the volume conduction model
for i=1:length(bnd)
  if ~isempty(map)
    edgecolor   = [color(i,1) color(i,2) color(i,3)];
    vertexcolor = [color(i,1) color(i,2) color(i,3)];
  end
  plot_mesh(bnd(i), 'faces', faces, 'vertices', vertices, 'edges', edges, 'facecolor', facecolor, 'facealpha', facealpha, 'vertexcolor', vertexcolor, 'edgecolor', edgecolor);
end
