function vol = pr_matrix2vol(mat3d, mat)
% returns (pseudo) vol struct for 3d matrix 
% FORMAT vol = pr_matrix2vol(mat3d,mat)
%
% Inputs
% mat3d   - 3D matrix
% mat     - optional 4x4 voxel -> world transformation
% 
% Outputs
% vol     - kind of SPM vol struct with matrix data added
%
% $Id: pr_matrix2vol.m,v 1.1 2005/04/20 15:05:00 matthewbrett Exp $
  
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
