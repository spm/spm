function [mat, dim] = spm_get_matdim(img, vx, bb)
% Voxel-to-world matrix and image dimensions from image or bbox and vox-dim
%
% FORMAT [mat, dim] = spm_get_matdim(img, vx, bb)
%
% img - filename of image to use as reference (defaults to SPM's TPM.nii)
% vx  - [1 x 3] vector of voxel dimensions (mm).
% bb  - [2 x 3] array of the min and max X, Y, and Z coordinates (mm),
%       i.e. bb = [minX minY minZ; maxX maxY maxZ].
%
% mat - [4 x 4] matrix mapping voxel coordinates to world (mm) coordinates
% dim - [1 x 3] vector of image dimensions (number of voxels)
%       (both as in output from spm_vol)
%
%       Note that the output mat will correspond to the same orientation
%       as SPM's canonical templates (transverse and vx(1) forced negative)
%       if either or both bb and vx are specified (finite), but otherwise
%       will keep the orientation of the reference image.
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: spm_get_matdim.m 5374 2013-03-29 17:26:24Z ged $

if nargin < 3, vx  = nan(1, 3); end
if nargin < 2, bb  = nan(2, 3); end
if nargin < 1, img = '';        end

% Use MNI space by default, based on tissue priors
try
    vol = spm_vol(img);
    if isempty(vol), error('Failed to read volume %s', img), end
catch
    vol = spm_vol(fullfile(spm('dir'), 'tpm', 'TPM.nii'));
end
mat = vol(1).mat;
dim = vol(1).dim;

valid_bb = all(isfinite(bb(:)));
valid_vx = all(isfinite(vx));
if ~valid_bb && ~valid_vx
    return
end

% User has specified one or both of bb or vx, over-ride appropriately
[BB VX] = spm_get_bbox(vol(1));
if ~valid_bb, bb = BB; end
if ~valid_vx, vx = VX; end

% Determine mat and dim from bb and vx, assuming "canonical" orientation
% (i.e. transverse with negative first voxel dimension)
vx  = [-1 1 1] .* abs(vx);
mn  = vx .* min(bb ./ repmat(vx, 2, 1)); % "first" voxel's mm coordinates
mx  = vx .* round(max(bb ./ repmat(vx, 2, 1))); % "last voxel's mm coords
% matrix that maps voxel [1 1 1] to mn
mat = spm_matrix([mn 0 0 0 vx]) * spm_matrix([-1 -1 -1]);
% dim such that mat * [dim 1]' == [mx 1]'
dim = mat \ [mx 1]';
dim = round(dim(1:3)');
