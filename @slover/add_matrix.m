function obj = add_matrix(obj, mat3d, mat, imgno)
% adds 3d matrix image vol to slice overlay 
% FORMAT obj = add_matrix(imgno, mat3d, mat)
%
% Inputs
% obj          - object
% mat3d        - 3D matrix to add as img
% mat          - optional 4x4 voxel->world translation
% imgno        - optional img no to add to (defaults to last in object)
% 
% Ouputs
% obj          - modified object
%
% $Id: add_matrix.m,v 1.1 2005/04/20 15:05:36 matthewbrett Exp $

if nargin < 2
  return
end
if nargin < 3
  mat = [];
end
if nargin < 4
  imgno = [];
end
if isempty(imgno)
  imgno = length(obj.img);
end
if ~isempty(mat3d)
  obj.img(imgno).vol = pr_matrix2vol(mat3d, mat);
end
