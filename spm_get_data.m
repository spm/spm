function [Y] = spm_get_data(V,XYZ)
% gets data from image files at specified locations
% FORMAT [Y] = spm_get_data(V,XYZ);
%
% V    - [1 x n] struct array of file handles (or filename matrix)
% XYZ  - [4 x m] location matrix (voxel)
%
% Y    - (n x m) double values
%
% see spm_sample_vol
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_get_data.m 112 2005-05-04 18:20:52Z john $


% ensure V is an array of handle structures
%---------------------------------------------------------------------------
if ~isstruct(V)
	V = spm_vol(V);
	try
		V = cat(2,V{:});
	end
end

% get data
%---------------------------------------------------------------
Y     = zeros(length(V),size(XYZ,2));
for i = 1:length(V)

	%-Load mask image within current mask & update mask
	%-------------------------------------------------------
	Y(i,:) = spm_sample_vol(V(i),XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
end

