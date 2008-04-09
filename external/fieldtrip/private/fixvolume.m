function [volume] = fixvolume(volume);

% FIXVOLUME ensures that a source reconstruction, MRI or statistical
% volume is described using a homogenous transformation matrix and that
% the volumes are reshaped to the right dimension.
%
% Furthermore, if the volume contains a seperate inside and outside
% vector with voxel indices, these are transformed into a single boolean 
% inside volume. 

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: fixvolume.m,v $
% Revision 1.2  2006/04/10 12:14:39  roboos
% use the new function private/fixinside instead of having a local snippet of code
%
% Revision 1.1  2006/02/27 16:55:55  roboos
% new function
%

% the volume should have a homogenous transformation matrix instead of xgrid/ygrid/zgrid
volume = grid2transform(volume);

% it should contain a logical ROI with the same dimensions as the volume itself
volume = fixinside(volume, 'logical');

% reshape each of the volumes that is found into a 3D array
param = parameterselection('all', volume);
for i=1:length(param)
  tmp    = getsubfield(volume, param{i});
  tmp    = reshape(tmp, volume.dim);
  volume = setsubfield(volume, param{i}, tmp);
end
