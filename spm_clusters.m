function A = spm_clusters(L,n)
% Return the cluster index for a point list
% FORMAT A = spm_clusters(L,n)
% L     - locations [x y x]' {in voxels} ([3 x m] matrix)
% n     - connectivity criterion (see spm_bwlabel) [Default: 18]
%
% A     - cluster index or region number ([1 x m] vector)
%__________________________________________________________________________
%
% spm_clusters characterises a point list of voxel values defined with
% their locations (L) in terms of edge, face and vertex connected
% subsets, returning a list of indices in A, such that the ith location
% belongs to cluster A(i) (using an 18 connectivity scheme).
%__________________________________________________________________________
% Copyright (C) 1994-2012 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson
% $Id: spm_clusters.m 4929 2012-09-17 14:21:01Z guillaume $


if isempty(L), A = []; return; end
if nargin < 2, n = 18; end

% Turn location list to binary 3D volume
%--------------------------------------------------------------------------
dim       = [max(L(1,:)) max(L(2,:)) max(L(3,:))];
vol       = zeros(dim(1),dim(2),dim(3));
indx      = sub2ind(dim,L(1,:)',L(2,:)',L(3,:)');
vol(indx) = 1;

% Label each cluster in 3D volume with its own label using an 18 
% connectivity criterion
%--------------------------------------------------------------------------
cci = spm_bwlabel(vol,n);

% Map back to list
%--------------------------------------------------------------------------
A = cci(indx);
A = A(:)';
