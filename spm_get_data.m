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
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_get_data.m 226 2005-09-12 15:00:59Z karl $


% ensure V is an array of handle structures
%--------------------------------------------------------------------------
if ~isstruct(V)
	V = spm_vol(V);
	try
		V = cat(2,V{:});
	end
end

% get data
%--------------------------------------------------------------------------
Y     = zeros(length(V),size(XYZ,2));
for i = 1:length(V)
    
        % check files exists, if not try pwd
    %----------------------------------------------------------------------
    if exist(V(i).fname,'file') ~=2
        [p,n,e]=fileparts(V(i).fname);
        V(i).fname = [n e];
    end

	%-Load mask image within current mask & update mask
	%----------------------------------------------------------------------
	Y(i,:) = spm_sample_vol(V(i),XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
end

