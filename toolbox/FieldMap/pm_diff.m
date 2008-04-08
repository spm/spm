function V = pm_diff(V,dir)
% Calculate derivative in one direction of 
% volume (matrix or memory mapped)
%
% FORMAT: V = pm_diff(V,dir)
%
% Input:
% V      : 3D matlab array, or filestruct returned
%          from spm_vol.
% dir    : Direction (1, 2 or 3 for x, y or z
%          respectively).
% Output:
% V      : 3D matlab array of derivatives.
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton
% $Id: pm_diff.m 1317 2008-04-08 16:16:38Z chloe $

if ischar(V) & exist(V) == 2
   V = spm_vol(V);
   dim = V.dim(1:3);
elseif isstruct(V) & isfield(V,'dim')
   dim = V.dim(1:3);
else
   dim=size(V);
end

hold =1;
[x,y,z]=ndgrid(1:dim(1),1:dim(2),1:dim(3));
[X,dX,dY,dZ] = spm_sample_vol(V,x,y,z,hold);

if dir==1
   V=reshape(dX,dim(1),dim(2),dim(3));
elseif dir==2
   V=reshape(dY,dim(1),dim(2),dim(3));
elseif dir==3
   V=reshape(dZ,dim(1),dim(2),dim(3));
end
