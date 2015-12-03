function vol = pr_matrix2vol(mat3d, mat)
% Return (pseudo) vol struct for 3d matrix
% FORMAT vol = pr_matrix2vol(mat3d,mat)
%
% Inputs
% mat3d   - 3D matrix
% mat     - optional 4x4 voxel -> world transformation
%
% Outputs
% vol     - kind of SPM vol struct with matrix data added
%__________________________________________________________________________

% Matthew Brett
% $Id: pr_matrix2vol.m 6623 2015-12-03 18:38:08Z guillaume $

if nargin < 1
    error('Need matrix to add to vol');
end
if nargin < 2
    mat = [];
end
if isempty(mat)
    mat = spm_matrix([]);
end
vol = [];
if ~isempty(mat3d)
    vol.imgdata = mat3d;
    vol.mat = mat;
    vol.dim = size(mat3d);
end
