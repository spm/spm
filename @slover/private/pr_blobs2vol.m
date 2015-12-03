function vol = pr_blobs2vol(xyz,vals,mat)
% Take XYZ matrix and values and return SPM matrix vol struct
% FORMAT vol = pr_blobs2vol(xyz,vals,mat)
%
% Inputs
% xyz      - 3xN X Y Z coordinate matrix (in voxels)
% vals     - 1xN values, one per coordinate
% mat      - 4x4 voxel->world space transformation
%
% Outputs
% vol      - vol struct, with matrix data 'imgdata' field
%__________________________________________________________________________

% Matthew Brett
% $Id: pr_blobs2vol.m 6623 2015-12-03 18:38:08Z guillaume $

vol = [];
if ~isempty(xyz),
    rcp      = round(xyz);
    vol.dim  = max(rcp,[],2)';
    off      = rcp(1,:) + vol.dim(1)*(rcp(2,:)-1+vol.dim(2)*(rcp(3,:)-1));
    vol.imgdata = zeros(vol.dim)+NaN;
    vol.imgdata(off) = vals;
    vol.imgdata      = reshape(vol.imgdata,vol.dim);
    vol.mat = mat;
end
