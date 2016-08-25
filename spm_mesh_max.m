function [N,Z,M,A,XYZ] = spm_mesh_max(X,L,G)
% Sizes, local maxima and locations of excursion sets on a surface mesh
% FORMAT [N,Z,M,A,XYZ] = spm_mesh_max(X,L,G)
% X        - a [nx1] array of stat values
% L        - a [nx1] array of locations {in vertices}
% G        - a patch structure
%
% N        - a [px1] size of connected components {in vertices}
% Z        - stat values of maxima
% M        - location of maxima {in vertices}
% A        - region number
% XYZ      - cell array of vertices locations
%__________________________________________________________________________
%
% See also: spm_max.m, spm_mesh_clusters.m and spm_mesh_get_lm.m
%__________________________________________________________________________
% Copyright (C) 2012-2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_max.m 6860 2016-08-25 12:00:10Z guillaume $


%-Get connected components
%--------------------------------------------------------------------------
LL     = NaN(size(G.vertices,1),1);
LL(L(1,:)) = X;
[C, N] = spm_mesh_clusters(G,LL);

%-Get local maxima
%--------------------------------------------------------------------------
M = spm_mesh_get_lm(G,LL);
Z = LL(M);
A = C(M);
M = [M;ones(2,size(M,2))];
N = N(A);

if nargout > 4
    XYZ = cell(1,max(A));
    for i=1:numel(XYZ)
        XYZ{i} = find(C==i)';
        XYZ{i} = [XYZ{i};ones(2,size(XYZ{i},2))];
    end
end
