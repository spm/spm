function [Y] = spm_get_data(V,XYZ,check)
% gets data from image files at specified locations
% FORMAT [Y] = spm_get_data(V,XYZ);
%
% V     - [1 x n] struct array of file handles (or filename matrix)
% XYZ   - [4 x m] or [3 x m]location matrix (voxel)
% check - check validity of input parameters [default: true]
%
% Y     - (n x m) double values
%
% see spm_sample_vol
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_get_data.m 3297 2009-07-29 17:20:19Z guillaume $

if nargin < 3, check = true; end

% ensure V is an array of handle structures
%--------------------------------------------------------------------------
if check && ~isstruct(V)
    V = spm_vol(V);
    try
        V = cat(2,V{:});
    end
end
 
% get data
%--------------------------------------------------------------------------
Y     = zeros(length(V),size(XYZ,2));
for i = 1:length(V)
 
    % check files exists and try pwd
    %----------------------------------------------------------------------
    if check && ~spm_existfile(V(i).fname)
        [p,n,e]     = fileparts(V(i).fname);
        V(i).fname  = [n,e];
    end
 
    %-Load data
    %----------------------------------------------------------------------
    Y(i,:) = spm_sample_vol(V(i),XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
 
end
