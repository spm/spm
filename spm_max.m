function [N,Z,M,A,XYZ] = spm_max(X,L)
% Sizes, maxima and locations of local excursion sets
% FORMAT [N Z M A XYZ] = spm_max(X,L)
% X     - values of 3-D field
% L     - locations [x y z]' {in voxels}
%
% N     - size of region {in voxels)
% Z     - Z values of maxima
% M     - location of maxima {in voxels}
% A     - region number
% XYZ   - cell array of voxel locations
%__________________________________________________________________________
%
% spm_max characterizes a point list of voxel values (X) and their
% locations (L) in terms of edge, face and vertex connected subsets,
% returning a maxima- orientated list:  The value of the ith maximum is
% Z(i) and its location is given by M(:,i). A(i) identifies the ith
% maximum with a region. Region A(i) contains N(i) voxels, whose
% coordinates are in a 3-by-N(i) array in XYZ{i}.
%
% See also: spm_bwlabel.m and spm_clusters.m
%__________________________________________________________________________
% Copyright (C) 2003-2011 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson
% $Id: spm_max.m 4384 2011-07-06 17:00:20Z guillaume $

if isempty(L)
    N = []; Z = []; M = []; A = []; XYZ = [];
    return;
end

%-Ensure that L contains exactly integers
%--------------------------------------------------------------------------
L          = round(L);

%-Turn location list to binary 3D volume
%--------------------------------------------------------------------------
dim        = [max(L(1,:)) max(L(2,:)) max(L(3,:))];
vol        = zeros(dim(1),dim(2),dim(3));
index      = sub2ind(dim,L(1,:)',L(2,:)',L(3,:)');
vol(index) = 1;

%-Label each cluster with its own label using an 18 connectivity criterion
% cci = connected components image volume
%--------------------------------------------------------------------------
[cci,num]  = spm_bwlabel(vol,18);

%-Get size (in no. of voxels) for each connected component
% ccs = connected component size
%--------------------------------------------------------------------------
ccs        = histc(cci(:),(0:num) + 0.5);
ccs        = ccs(1:end-1);

%-Get indices into L for voxels that are indeed local maxima (using an 18 
% neighbour criterion)
%--------------------------------------------------------------------------
vol(index) = X;
Lindex     = spm_get_lm(vol,L);

M          = L(:,Lindex);
Z          = X(Lindex);
mindex     = sub2ind(dim,L(1,Lindex)',L(2,Lindex)',L(3,Lindex)');
A          = cci(mindex);
N          = ccs(A);

%-Cell array of XYZ locations of voxels in each cluster
%--------------------------------------------------------------------------
if nargout > 4
    xyz(:,index) = sparse(L);
    cci   = sparse(cci(:));
    XYZ = cell(1, max(A));
    for i = 1:max(A)
        XYZ{i} = full(xyz(:,cci == i));
    end
end
