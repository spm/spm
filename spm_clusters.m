function A = spm_clusters(L)
% Return the cluster index for a point list
% FORMAT [A] = spm_clusters(L)
% L     - locations [x y x]' {in voxels} ([3 x m] matrix)
%
% A     - cluster index or region number ([1 x m] vector)
%__________________________________________________________________________
%
% spm_clusters characterizes a point list of voxel values defined with
% their locations (L) in terms of edge, face and vertex connected
% subsets, returning a list of indices in A, such that the ith location
% belongs to cluster A(i) (using an 18 connectivity scheme).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson
% $Id: spm_clusters.m 3822 2010-04-16 18:43:08Z karl $


if isempty(L), A = []; return; end


% Turn location list to binary 3D volume
%--------------------------------------------------------------------------
dim       = [max(L(1,:)) max(L(2,:)) max(L(3,:))];
vol       = zeros(dim(1),dim(2),dim(3));
indx      = sub2ind(dim,L(1,:)',L(2,:)',L(3,:)');
vol(indx) = 1;


% Label each cluster in 3D volume with its own label using an 18 
% connectivity criterion
%--------------------------------------------------------------------------
[cci,num] = spm_bwlabel(vol,18);


% Map back to list
%--------------------------------------------------------------------------
A = cci(indx);
A = A(:)';
