function [C, N] = spm_mesh_clusters(M,T)
% Label connected components of a surface mesh texture
% FORMAT [C, N] = spm_mesh_clusters(M,T)
% M        - a [mx3] faces array or a patch structure
% T        - a [nx1] texture vector (using NaN)
%
% C        - a [nx1] vector of cluster index
% N        - a [px1] size of connected components {in vertices}
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_clusters.m 3989 2010-07-12 18:39:16Z guillaume $


%-Input parameters
%--------------------------------------------------------------------------
if isnumeric(M)
    F   = M;
else
    F   = M.faces;
end

%-Compute a reduced mesh corresponding to the texture
%--------------------------------------------------------------------------
B       = NaN(size(T));
i       = find(~isnan(T));
B(i)    = 1:length(i);
F       = B(F);
F(any(isnan(F),2),:) = [];

%-Compute its adjacency matrix
%--------------------------------------------------------------------------
A       = spm_mesh_adjacency(F);
A       = A + speye(size(A));

%-Extract connected components using Dulmage-Mendelsohn decomposition
%--------------------------------------------------------------------------
[p,q,r] = dmperm(A);

CC      = T(i);
for j=1:length(r)-1
    CC(p(r(j):r(j+1)-1)) = j;
end

C       = NaN(numel(T),1);
C(i)    = CC;

%-Sort connected component labels according to their size
%--------------------------------------------------------------------------
N       = diff(r);
[N,ni]  = sort(N(:),1,'descend');
[ni,ni] = sort(ni);
C(i)    = ni(C(i));
