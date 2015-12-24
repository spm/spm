function D = pm_diff(V,dir)
% Calculate derivative in one direction of volume (matrix or memory mapped)
% FORMAT D = pm_diff(V,dir)
% V    - 3D array, or filestruct returned from spm_vol
% dir  - direction (1, 2 or 3 for x, y or z respectively)
%
% D    - 3D array of derivatives
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton
% $Id: pm_diff.m 6656 2015-12-24 16:49:52Z guillaume $

if ischar(V)
   V   = spm_vol(V);
end
if isstruct(V)
   dim = V.dim;
else
   dim = size(V);
end
dim    = [dim 1 1 1];

hold   = 1;
[x,y,z]      = ndgrid(1:dim(1),1:dim(2),1:dim(3));
[X,dX,dY,dZ] = spm_sample_vol(V,x,y,z,hold);

switch dir
	case 1
        D    = reshape(dX,dim(1),dim(2),dim(3));
    case 2
        D    = reshape(dY,dim(1),dim(2),dim(3));
    case 3
        D    = reshape(dZ,dim(1),dim(2),dim(3));
    otherwise
        error('Unknown direction.');
end
