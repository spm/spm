function obj = add_blobs(obj, xyz, vals, mat, imgno)
% Add SPM blobs to img no 'imgno', as specified in
% FORMAT obj = add_blobs(obj, xyz, vals, mat, imgno)
%
% Inputs
% XYZ   - 3xN voxel coordinates of N blob values
% vals  - N blob intensity values
% mat   - 4x4 matrix specifying voxels -> mm
% imgno - slice overlay img number to add to (defaults last in object)
%
% Outputs
% obj   - modified object
%__________________________________________________________________________

% Matthew Brett
% $Id: add_blobs.m 6623 2015-12-03 18:38:08Z guillaume $

if nargin < 4
    error('Need all of object, xyz, vals, mat');
end
if nargin < 5
    imgno = [];
end
if isempty(imgno)
    imgno = length(obj.img);
end
if ~isempty(xyz)
    obj.img(imgno).vol = pr_blobs2vol(xyz,vals,mat);
end
