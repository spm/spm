function A = spm_clusters(L)
% Returns the cluster index for a point list
% FORMAT [A] = spm_clusters(L)
% L     - locations [x y x]' {in voxels}
%
% A     - cluster index or region number
%_______________________________________________________________________
%
% spm_clusters characterizes a point list of voxel values defined with
% their locations (L) in terms of edge, face and vertex connected
% subsets, returning a list of indices in A, such that the ith location
% belongs to cluster A(i) (using an 18 connectivity scheme).
%
% The programming interface of this routine is modeled on the "old"
% spm_clusters to avoid the need to recode any other parts. The "old"
% version had the unfourtunate tendency to crash SPM (and Matlab) 
% whenever invoked with a "too" long list of locations. This new
% version has the same functionality, but because it is not recursive 
% it will (hopefully) not crash.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jesper Andersson
% $Id: spm_clusters.m 112 2005-05-04 18:20:52Z john $


%
% Turn location list to binary 3D volume.
%
dim = [max(L(1,:)) max(L(2,:)) max(L(3,:))];
vol = zeros(dim(1),dim(2),dim(3));
indx = sub2ind(dim,L(1,:)',L(2,:)',L(3,:)');
vol(indx) = 1;

%
% Label each cluster in 3D volume with it's 
% own little label using an 18 connectivity
% criterion (without crashing ;-)).
%

[cci,num] = spm_bwlabel(vol,18);

%
% Map back to list.
%

A = cci(indx');

return

