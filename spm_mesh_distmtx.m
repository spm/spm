function D = spm_mesh_distmtx(M,order)
% Compute the distance matrix of a triangle mesh
% FORMAT D = spm_mesh_distmtx(M,order)
% M        - patch structure
% order    - 0: adjacency matrix
%            1: first order distance matrix [default]
%            2: second order distance matrix
%
% D        - distance matrix
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_distmtx.m 4035 2010-08-05 18:54:32Z guillaume $

if nargin < 2, order = 1; end

if order > 2, error('High order distance matrix not handled.'); end

%-Adjacency matrix
%--------------------------------------------------------------------------
if order == 0
    D = spm_mesh_adjacency(M);
    return;
end

%-First order distance matrix
%--------------------------------------------------------------------------
if isstruct(M) && ~isa(M.faces,'double')
    M.faces = double(M.faces);
elseif isa(M,'gifti')
    M = export(M,'patch');
    M.faces = double(M.faces);
end
d = M.vertices(M.faces(:,[1 2 3]),:) - M.vertices(M.faces(:,[2 3 1]),:);

D = sparse([M.faces(:,1); M.faces(:,2); M.faces(:,3)], ...
           [M.faces(:,2); M.faces(:,3); M.faces(:,1)], sqrt(sum(d.^2, 2)));

D = (D + D')/2;

if order == 1, return; end

%-Second order distance matrix
%--------------------------------------------------------------------------
D2 = D;
for i=1:size(M.vertices,1)
    a      = find(D2(i,:));
    [b,c]  = find(D2(a,:));
    D(i,c) = D(i,a(b)) + diag(D(a(b),c))';
    D(c,i) = D(i,c)';
end
