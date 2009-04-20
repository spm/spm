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
%     'faces'         ['yes', 'no', 1, 0, 'true', 'false']
%     'facecolor'     ['brain', 'cortex', 'skin', 'black', 'red', 'r', ..., [0.5 1 0], ...]
%     'faceindex'     ['yes', 'no', 1, 0, 'true', 'false']
%     'vertices'      ['yes', 'no', 1, 0, 'true', 'false']
%     'vertexcolor'   ['brain', 'cortex', 'skin', 'black', 'red', 'r', ..., [0.5 1 0], ...]
%     'vertexindex'   ['yes', 'no', 1, 0, 'true', 'false']
%     'edges'         ['yes', 'no', 1, 0, 'true', 'false']
%     'edgecolor'     ['brain', 'cortex', 'skin', 'black', 'red', 'r', ..., [0.5 1 0], ...] 
%     'colormap'      ['gray', 'hot', 'cool', 'copper', 'spring', 'summer', ...]
%
% Example
%   vol.r = [1 5 10];
%   vol.o = [0 0 4];
%   figure, plot_vol(vol)
%
% Copyright (C) 2009, Cristiano Micheli 
%
% $Log: plot_vol.m,v $
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
faces       = keyval('faces',      varargin);
facecolor   = keyval('facecolor',  varargin);
faceindex   = keyval('faceindex',  varargin);
vertices    = keyval('vertices',   varargin);  if isempty(vertices),vertices=1;end
vertexcolor = keyval('vertexcolor',  varargin);
vertexindex = keyval('vertexindex',  varargin);
vertexsize  = keyval('vertexsize',  varargin); if isempty(vertexsize),vertexsize=10;end
edges       = keyval('edges',      varargin);  if isempty(edges),edges=1;end
edgecolor   = keyval('edgecolor',      varargin);  
map         = keyval('colormap',      varargin); 

% we will probably need a sphere, so let's prepare one
[pnt, tri] = icosahedron642;

if ~isempty(map)
  try
    cmap=colormap(map);close(gcf);
  catch
    error('Colormap does not exist')
  end
  % set the color of the sphere:
  color = [(linspace(cmap(1,1),cmap(64,1),length(vol.r)))' (linspace(cmap(1,2),cmap(64,2),length(vol.r)))' ...
              (linspace(cmap(1,3),cmap(64,3),length(vol.r)))'];
end

switch voltype(vol)   
  case {'singlesphere' 'concentric'}

    vol.r = sort(vol.r);
    for i=1:length(vol.r)
      bnd = [];
      bnd.pnt = pnt*vol.r(i);
      bnd.pnt(:,1) = bnd.pnt(:,1) + vol.o(1);
      bnd.pnt(:,2) = bnd.pnt(:,2) + vol.o(2);
      bnd.pnt(:,3) = bnd.pnt(:,3) + vol.o(3);

      bnd.tri = tri;

      if ~isempty(map)
        plot_mesh(bnd,'edgecolor', [color(i,1) color(i,2) color(i,3)], ...
                        'vertexcolor', [color(i,1) color(i,2) color(i,3)]);
      else
        plot_mesh(bnd);
      end
    end
    
  case 'multisphere'
    for i=1:length(vol.label)
      bnd = [];
      bnd.pnt = pnt*vol.r(i);
      bnd.pnt(:,1) = bnd.pnt(:,1) + vol.o(i,1);
      bnd.pnt(:,2) = bnd.pnt(:,2) + vol.o(i,2);
      bnd.pnt(:,3) = bnd.pnt(:,3) + vol.o(i,3);

      bnd.tri = tri;
      
      if ~isempty(map)
        plot_mesh(bnd,'edgecolor', [color(i,1) color(i,2) color(i,3)], ...
                        'vertexcolor', [color(i,1) color(i,2) color(i,3)]);
      else
        plot_mesh(bnd);
      end
    end    

  case {'bem', 'dipoli', 'asa', 'avo', 'bemcp', 'nolte'}
    for i=1:length(vol.bnd)
      if ~isempty(map)
        plot_mesh(vol.bnd(i),'edgecolor', [color(i,1) color(i,2) color(i,3)], ...
                        'vertexcolor', [color(i,1) color(i,2) color(i,3)]);
      else
        plot_mesh(vol.bnd(i));
      end      
    end
   
  otherwise
    error('unsupported voltype')
end

