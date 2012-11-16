function [C, N] = spm_mesh_clusters(M,T)
% Label connected components of surface mesh data
% FORMAT [C, N] = spm_mesh_clusters(M,T)
% M        - a [mx3] faces array or a patch structure
% T        - a [nx1] data vector (using NaNs or logicals), n = #vertices
%
% C        - a [nx1] vector of cluster indices
% N        - a [px1] size of connected components {in vertices}
%__________________________________________________________________________
% Copyright (C) 2010-2012 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_clusters.m 5065 2012-11-16 20:00:21Z guillaume $


%-Input parameters
%--------------------------------------------------------------------------
if ~islogical(T)
    T   = ~isnan(T);
end

%-Compute the (reduced) adjacency matrix
%--------------------------------------------------------------------------
A       = spm_mesh_adjacency(M);
A       = A + speye(size(A));
A(~T,:) = [];
A(:,~T) = [];

%-And perform Dulmage-Mendelsohn decomposition to find connected components
%--------------------------------------------------------------------------
[p,q,r] = dmperm(A);
N       = diff(r);
CC      = zeros(size(A,1),1);
for i=1:length(r)-1
    CC(p(r(i):r(i+1)-1)) = i;
end
C       = NaN(numel(T),1);
C(T)    = CC;

%-Sort connected component labels according to their size
%--------------------------------------------------------------------------
[N,ni]  = sort(N(:), 1, 'descend');
[ni,ni] = sort(ni);
C(T)    = ni(C(T));
