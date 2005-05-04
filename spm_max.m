function [N,Z,M,A] = spm_max(X,L)
% Sizes, maxima and locations of local excursion sets
% FORMAT [N Z M A] = spm_max(X,L)
% X     - values of 3-D field
% L     - locations [x y x]' {in voxels}
% N     - size of region {in voxels)
% Z     - Z values of maxima
% M     - location of maxima {in voxels}
% A     - region number
%_______________________________________________________________________
%
% spm_max characterizes a point list of voxel values (X) and their
% locations (L) in terms of edge, face and vertex connected subsets,
% returning a maxima- orientated list:  The value of the ith maximum is
% Z(i) and its location is given by M(:,i). A(i) identifies the ith
% maximum with a region. Region A(i) contains N(i) voxels.
%
% See also: spm_bwlabel.c and spm_clusters.m
%
% Rewrite of "old" spm_max (that used a recursive algorithm) for the
% connected-component labelling that is part of determining clusters.
% The programming interface is identical to the old to avoid having
% to recode any routines calling spm_max or spm_clusters.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jesper Andersson
% $Id: spm_max.m 112 2005-05-04 18:20:52Z john $


% Ensure that L contains exactly integers
L = round(L);

%
% Turn location list to binary 3D volume.
%
dim = [max(L(1,:)) max(L(2,:)) max(L(3,:))];
vol = zeros(dim(1),dim(2),dim(3));
index = sub2ind(dim,L(1,:)',L(2,:)',L(3,:)');
vol(index) = 1;

%
% Label each cluster in 3D volume with it's 
% own little label using an 18 connectivity
% criterion (without crashing ;-)).
%

% cci = connected components image volume.

[cci,num] = spm_bwlabel(vol,18);

%
% Get size (in no. of voxels) for each
% connected component.
%

% ccs = connected component size

ccs = histc(cci(:),[0:max(cci(:))]+0.5);
ccs = ccs(1:end-1);

%
% Get indicies into L for voxels that are
% indeed local maxima (using an 18 neighbour
% criterion).
%

vol(index) = X;
Lindex = spm_get_lm(vol,L); % That was soo not intended.

M = L(:,Lindex);
Z = X(Lindex);
mindex = sub2ind(dim,L(1,Lindex)',L(2,Lindex)',L(3,Lindex)');
A = cci(mindex);
N = ccs(A);

return

