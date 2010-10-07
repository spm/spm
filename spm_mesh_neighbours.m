function [N,D] = spm_mesh_neighbours(M,order)
% Return first-order neighbours of all vertices of a surface mesh
% FORMAT N = spm_mesh_neighbours(M,order)
% M        - a patch structure or an adjacency matrix
% order    - ordinal or euclidean distance for 1st order neighbours {[0],1}
%
% N        - a [nxp] neighbours array (n = #vertices, p = # max neighbours)
% D        - a [nxp] distance array to neighbours
% N & D contain 0 when number of neighbours is smaller than p.
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_neighbours.m 4081 2010-10-07 14:04:44Z guillaume $

if nargin < 2, order = 0; end

%-Compute adjacency matrix
%--------------------------------------------------------------------------
if issparse(M) && size(M,1) == size(M,2)
    A = M;
else
    A = spm_mesh_distmtx(M,order);
end

%-Build neighbours array
%--------------------------------------------------------------------------
if nargout > 1
    [N, D] = spm_mesh_utils('neighbours',A);
else
    N = spm_mesh_utils('neighbours',A);
end

% [i, j, s] = find(A);
% p = max(sum(A~=0,2));
% N = zeros(size(A,1),p);
% if nargout > 1, D = zeros(size(N)); end
% l = 1;
% for k=1:length(j)
%     N(j(k),l) = i(k);
%     if nargout > 1, D(j(k),l) = s(k); end
%     try, if j(k+1) ~= j(k), l = 1; else l = l + 1; end; end
% end
