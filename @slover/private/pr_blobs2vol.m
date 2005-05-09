function vol = pr_blobs2vol(xyz,vals,mat)
% takes XYZ matrix and values, returns SPM matrix vol struct
% FORMAT vol = pr_blobs2vol(xyz,vals,mat)
% 
% Inputs 
% xyz      - 3xN X Y Z coordinate matrix (in voxels)
% vals     - 1xN values, one per coordinate
% mat      - 4x4 voxel->world space transformation
% 
% Outputs
% vol      - vol struct, with matrix data 'imgdata' field
% 
% $Id: pr_blobs2vol.m,v 1.1 2005/04/20 15:05:00 matthewbrett Exp $
  
if nargin < 3
  error('Need XYZ, vals and mat');
end

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
return

