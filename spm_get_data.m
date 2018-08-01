function [Y] = spm_get_data(V,XYZ,check)
% Get data from image files at specified locations
% FORMAT [Y] = spm_get_data(V,XYZ,check)
%
% V          - [1 x n] struct array of file handles (or filename matrix)
% XYZ        - [3 x m] or [4 x m] location matrix {voxel}
% check      - check validity of input parameters [default: true]
%
% Y          - [n x m] double values
%
% See also spm_sample_vol
%__________________________________________________________________________
% Copyright (C) 2002-2018 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_get_data.m 7385 2018-08-01 15:58:24Z guillaume $


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
Y     = zeros(numel(V),size(XYZ,2));
for i = 1:numel(V)
 
    %-Check files exists and try current directory otherwise
    %----------------------------------------------------------------------
    fname           = V(i).fname;
    if check && ~spm_existfile(fname)
        [p,n,e]     = fileparts(fname);
        V(i).fname  = [n,e];
    end
 
    %-Load data
    %----------------------------------------------------------------------
    try
        Y(i,:) = spm_sample_vol(V(i),XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
    catch
       fprintf('Could not access file "%s".\n',fname);
       rethrow(lasterror);
    end
 
end
