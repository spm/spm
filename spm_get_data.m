function [Y] = spm_get_data(V,XYZ,check)
% Get data from image files at specified locations
% FORMAT [Y] = spm_get_data(V,XYZ,check)
%
% V     - [1 x n] struct array of file handles (or filename matrix)
% XYZ   - [3 x m] or [4 x m] location matrix (voxel)
% check - check validity of input parameters [default: true]
%
% Y     - (n x m) double values
%
% See also spm_sample_vol
%__________________________________________________________________________
% Copyright (C) 2002-2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_get_data.m 5160 2012-12-21 16:58:38Z guillaume $

if nargin < 3, check = true; end

%-Ensure V is an array of handle structures
%--------------------------------------------------------------------------
if check && ~isstruct(V)
    V = spm_vol(V);
    try
        V = cat(2,V{:});
    end
end
 
%-Get data
%--------------------------------------------------------------------------
Y     = zeros(length(V),size(XYZ,2));
for i = 1:length(V)
 
    %-Check files exists and try current directory otherwise
    %----------------------------------------------------------------------
    if check && ~spm_existfile(V(i).fname)
        [p,n,e]     = fileparts(V(i).fname);
        V(i).fname  = [n,e];
    end
 
    %-Load data
    %----------------------------------------------------------------------
    try
        Y(i,:) = spm_sample_vol(V(i),XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
    catch
       fprintf('Error with file %s\n',V(i).fname);
       rethrow(lasterror);
    end
 
end
