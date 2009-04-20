function [hs, hc, contour] = plot_mesh(bnd, varargin)

% PLOT_MESH visualizes the information of a mesh contained in the first
% argument bnd. The boundary argument (bnd) contains typically 2 fields
% called .pnt and .tri referring to vertices and triangulation of a mesh.
%
% Use as
%   hs = plot_mesh(bnd, varargin)
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
%
% Example
%   [pnt, tri] = icosahedron162;
%   bnd.pnt = pnt;
%   bnd.tri = tri;
%   plot_mesh(bnd, 'faces', 'yes', 'vertices', 'yes', 'edges', 'no', 'facecolor', 'skin')

% Copyright (C) 2009, Cristiano Micheli
%
% $Log: plot_mesh.m,v $
% Revision 1.10  2009/04/17 13:43:33  crimic
% updated help
%
% Revision 1.9  2009/04/17 13:38:23  crimic
% added default options for varargin, added edgecolor argument
%
% Revision 1.8  2009/04/09 09:36:50  crimic
% *** empty log message ***
%
% Revision 1.7  2009/04/09 09:35:53  crimic
% integrated help
%
% Revision 1.6  2009/04/09 09:20:33  crimic
% clean-up of non used options (val, contour and surface), inserted istrue check
%
% Revision 1.5  2009/04/08 17:09:40  crimic
% added help and modified input argument structure
%
% Revision 1.4  2009/04/08 16:08:50  crimic
% indented with 2 spaces tab, added log signature and copyright
%

% FIXME: introduce option for color coding (see sourceplot)


if ~isfield(bnd,'pnt') && ~isfield(bnd,'tri')
  error('First argument must be a boundary structure (points, triangulation)')
end

% get the optional input arguments
faces       = keyval('faces',      varargin);
facecolor   = keyval('facecolor',  varargin);
faceindex   = keyval('faceindex',  varargin);
vertices    = keyval('vertices',   varargin);  if isempty(vertices),vertices=1;end
vertexcolor = keyval('vertexcolor',  varargin);
vertexindex = keyval('vertexindex',  varargin);
vertexsize  = keyval('vertexsize',  varargin); if isempty(vertexsize),vertexsize=10;end
edges       = keyval('edges',      varargin);  if isempty(edges),edges=1;end
edgecolor   = keyval('edgecolor',      varargin);  if isempty(edgecolor),edgecolor='k';end

% start with empty return values
hs      = [];
skin   = [255 213 119]/255;
skull  = [140  85  85]/255;
brain  = [202 100 100]/255;
cortex = [255 213 119]/255;

% new colors management
if strcmpi(vertexcolor,'skin') || strcmpi(vertexcolor,'brain') || strcmpi(vertexcolor,'cortex')
  vertexcolor = eval(vertexcolor);
end
if strcmpi(facecolor,'skin') || strcmpi(facecolor,'brain') || strcmpi(facecolor,'cortex')
  facecolor = eval(facecolor);
end

% everything is added to the current figure
holdflag = ishold;
hold on

pnt = bnd.pnt;
tri = bnd.tri;

if istrue(faces)
  hs = patch('Vertices', pnt, 'Faces', tri);
  set(hs, 'FaceColor', 'white');
  set(hs, 'EdgeColor', 'none');
  if ~isempty(facecolor)
    try
      set(hs, 'FaceColor', facecolor);
    catch
      error('Unknown color')
    end
  end
  if istrue(faceindex)
    % plot the triangle indices (numbers) at each face
    for face_indx=1:size(tri,1)
      str = sprintf('%d', face_indx);
      tri_x = (pnt(tri(face_indx,1), 1) +  pnt(tri(face_indx,2), 1) +  pnt(tri(face_indx,3), 1))/3;
      tri_y = (pnt(tri(face_indx,1), 2) +  pnt(tri(face_indx,2), 2) +  pnt(tri(face_indx,3), 2))/3;
      tri_z = (pnt(tri(face_indx,1), 3) +  pnt(tri(face_indx,2), 3) +  pnt(tri(face_indx,3), 3))/3;
      h   = text(tri_x, tri_y, tri_z, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
      hs  = [hs; h];
    end
  end
end

if istrue(vertices)
  if size(pnt, 2)==2
    hs = plot(pnt(:,1), pnt(:,2), 'k.');
  else
    hs = plot3(pnt(:,1), pnt(:,2), pnt(:,3), 'k.');
  end
  if ~isempty(vertexcolor)
    try
      set(hs, 'Marker','.','MarkerEdgeColor', vertexcolor,'MarkerSize', vertexsize);
    catch
      error('Unknown color')
    end
  end
  if istrue(vertexindex)
    % plot the vertex indices (numbers) at each node
    for node_indx=1:size(pnt,1)
      str = sprintf('%d', node_indx);
      if size(pnt, 2)==2
        h   = text(pnt(node_indx, 1), pnt(node_indx, 2), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
      else
        h   = text(pnt(node_indx, 1), pnt(node_indx, 2), pnt(node_indx, 3), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
      end
      hs  = [hs; h];
    end
  end
end

if istrue(edges)
  % plot the edges of the 2D or 3D triangulation
  hs = patch('Vertices', pnt, 'Faces', tri);
  set(hs, 'FaceColor', 'none');
  set(hs, 'EdgeColor', edgecolor);
end

axis off
axis vis3d
axis equal

if nargout==0
  clear hs
end
if ~holdflag
  hold off
end

