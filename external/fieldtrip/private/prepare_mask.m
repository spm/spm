function [mask] = prepare_mask(layout, varargin)
% PREPARE_MASK  calculate anatomical mask only once, based only on layout
% Use as:
%   [mask] = prepare_mask(layout, 'key', 'val');
% with optional parameter:
% 'gridscale'    scale of the interpolated grid (default: 67)
%
% See also COMPONENTBROWSER, DATABROWSER

% Copyright (C) 2009
%
% $Log: prepare_mask.m,v $
% Revision 1.1  2009/07/15 08:40:54  giopia
% from topoplot, it relies on inside_contour
%

% check the input
if ~isfield(layout, 'mask')
  warning('prepare_mask:nomask', 'Mask is not present in layout. Mask not created');
  mask = [];
  return
end

[gridscale] = keyval('gridscale', varargin);
if isempty(gridscale);  gridscale = 67;  end

% find limits for interpolation:
xmin = +inf;
xmax = -inf;
ymin = +inf;
ymax = -inf;

for i=1:length(layout.mask)
  xmin = min([xmin; layout.mask{i}(:,1)]);
  xmax = max([xmax; layout.mask{i}(:,1)]);
  ymin = min([ymin; layout.mask{i}(:,2)]);
  ymax = max([ymax; layout.mask{i}(:,2)]);
end

xi = linspace(xmin, xmax, gridscale);   % x-axis for interpolation (row vector)
yi = linspace(ymin, ymax, gridscale);   % y-axis for interpolation (row vector)
Xi =  ones(gridscale,1)*xi;
Yi = (ones(gridscale,1)*yi)';

% apply anatomical mask to the data, i.e. that determines that the interpolated data outside the circle is not displayed
mask = false(gridscale);
for i=1:length(layout.mask)
  layout.mask{i}(end+1,:) = layout.mask{i}(1,:); % force them to be closed
  mask(inside_contour([Xi(:) Yi(:)], layout.mask{i})) = true;
end
