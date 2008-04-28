function [shape] = transform_headshape(transform, shape)

% TRANSFORM_HEADSHAPE applies a homogenous coordinate transformation to a
% structure with headshape and fiducial information.
%
% Use as
%   shape = transform_headshape(transform, shape)
%
% See also READ_HEADSHAPE

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: transform_headshape.m,v $
% Revision 1.2  2008/04/15 15:33:58  roboos
% fixed small bug (thanks to Vladimir)
%
% Revision 1.1  2008/04/11 16:14:03  roboos
% first implementation, simple helper function for spm integration and symmetry with vol and sens
%

if any(transform(4,:) ~= [0 0 0 1])
  error('invalid transformation matrix');
end

if isfield(shape, 'pnt') && ~isempty(shape.pnt)
  shape.pnt = apply(transform, shape.pnt);
end

if isfield(shape, 'fid') && isfield(shape.fid, 'pnt')
  shape.fid.pnt = apply(transform, shape.fid.pnt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that applies the homogenous transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new] = apply(transform, old)
old(:,4) = 1;
new = old * transform';
new = new(:,1:3);
