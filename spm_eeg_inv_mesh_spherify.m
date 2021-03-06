function [pnt, tri] = spm_eeg_inv_mesh_spherify(pnt, tri, varargin)
% Takes a cortical mesh and scales it so that it fits into a
% unit sphere.
%
% This function determines the points of the original mesh that support a
% convex hull and determines the radius of those points. Subsequently the
% radius of the support points is interpolated onto all vertices of the
% original mesh, and the vertices of the original mesh are scaled by
% dividing them by this interpolated radius.
%
% Use as
%   [pnt, tri] = spm_eeg_inv_mesh_spherify(pnt, tri, ...)
%
% Optional arguments should come as key-value pairs and may include
%   shift  = 'no', mean', 'range'
%   smooth = number (default = 20)
%__________________________________________________________________________

% Robert Oostenveld
% Copyright (C) 2009-2022 Wellcome Centre for Human Neuroimaging


% graphical display for debugging
fb = false;

shift  = keyval('shift', varargin);
smooth = keyval('smooth', varargin);

% set the concentration factor
if ~isempty(smooth)
  k = smooth;
else
  k = 100;
end

% the following code is for debugging
if fb
  figure
  [sphere_pnt, sphere_tri] = icosahedron162;
  y = vonmisesfischer(5, [0 0 1], sphere_pnt);
  triplot(sphere_pnt, sphere_tri, y);
end

npnt = size(pnt, 1);
ntri = size(tri, 1);

switch shift
  case 'mean'
    pnt(:,1) = pnt(:,1) - mean(pnt(:,1));
    pnt(:,2) = pnt(:,2) - mean(pnt(:,2));
    pnt(:,3) = pnt(:,3) - mean(pnt(:,3));
  case 'range'
    minx = min(pnt(:,1));
    miny = min(pnt(:,2));
    minz = min(pnt(:,3));
    maxx = max(pnt(:,1));
    maxy = max(pnt(:,2));
    maxz = max(pnt(:,3));
    pnt(:,1) = pnt(:,1) - mean([minx maxx]);
    pnt(:,2) = pnt(:,2) - mean([miny maxy]);
    pnt(:,3) = pnt(:,3) - mean([minz maxz]);
  otherwise
    % do nothing
end

% determine the convex hull, especially to determine the support points
tric = convhulln(pnt);
sel  = unique(tric(:));

% create a triangulation for only the support points
support_pnt = pnt(sel,:);
support_tri = convhulln(support_pnt);

if fb
  figure
  triplot(support_pnt, support_tri, [], 'faces_skin');
  triplot(pnt, tri, [], 'faces_skin');
  alpha 0.5
end

% determine the radius and thereby scaling factor for the support points
support_scale = zeros(length(sel),1);
for i=1:length(sel)
  support_scale(i) = norm(support_pnt(i,:));
end

% interpolate the scaling factor for the support points to all points
scale = zeros(npnt,1);
for i=1:npnt
  u = pnt(i,:);
  y = vonmisesfischer(k, u, support_pnt);
  y = y ./ sum(y);
  scale(i) = y' * support_scale;
end

% apply the interpolated scaling to all points
pnt(:,1) = pnt(:,1) ./ scale;
pnt(:,2) = pnt(:,2) ./ scale;
pnt(:,3) = pnt(:,3) ./ scale;

% downscale the points further to make sure that nothing sticks out
n = zeros(npnt,1);
for i = 1:npnt
    n(i) = norm(pnt(i, :));
end
mscale = (1-eps) / max(n);
pnt = mscale * pnt;

if fb
  figure
  [sphere_pnt, sphere_tri] = icosahedron162;
  triplot(sphere_pnt, sphere_tri, [], 'faces_skin');
  triplot(pnt, tri, [], 'faces_skin');
  alpha 0.5
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VONMISESFISCHER probability distribution
%
% Use as
%   y = vonmisesfischer(k, u, x)
% where
%   k = concentration parameter
%   u = mean direction
%   x = direction of the points on the sphere
%
% The von Mises?Fisher distribution is a probability distribution on the
% (p?1) dimensional sphere in Rp. If p=2 the distribution reduces to the
% von Mises distribution on the circle. The distribution belongs to the
% field of directional statistics.
%
% This implementation is based on
% http://en.wikipedia.org/wiki/Von_Mises-Fisher_distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = vonmisesfischer(k, u, x)

% the data describes N points in P-dimensional space
[n, p] = size(x);

% ensure that the direction vectors are unit length
u = u ./ norm(u);
for i=1:n
  x(i,:) = x(i,:) ./ norm(x(i,:));
end

% FIXME this normalisation is wrong
% but it is not yet needed, so the problem is acceptable for now
% Cpk = (k^((p/2)-1)) ./ ( (2*pi)^(p/2) * besseli(p/2-1, k));
Cpk = 1;

y = exp(k * u * x') ./ Cpk;
y = y(:);


function [val] = keyval(key, varargin)

% KEYVAL returns the value that corresponds to the requested key in a
% key-value pair list of variable input arguments
%
% Use as
%   [val] = keyval(key, varargin)
%
% See also VARARGIN

% Copyright (C) 2005-2007, Robert Oostenveld


if length(varargin)==1 && iscell(varargin{1})
  varargin = varargin{1};
end

if mod(length(varargin),2)
  error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
end

keys = varargin(1:2:end);
vals = varargin(2:2:end);

hit = find(strcmp(key, keys));
if length(hit)==0
  % the requested key was not found
  val = [];
elseif length(hit)==1  
  % the requested key was  found
  val = vals{hit};
else
  error('multiple input arguments with the same name');
end
