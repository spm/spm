function [Mdist] = spm_eeg_inv_meshdist(vert,face,order)
% Efficient computation of the 2nd order distance matrix of a triangulated
% irregular mesh, based on the cortical neighbourhoud (geodesic distance)
% vert  - locations of vertices
% face  - indices of faces
% order - 0; Adjacency matrix
%         1: 1st order distance matrix
%         2: 2nd order distance matrix
% Inspired by function mesh_laplacian.m by Darren Weber
% from the bioelectromagnetism matlab toolbox
% see http://eeg.sourceforge.net/

%==========================================================================

% default order is 2nd
%--------------------------------------------------------------------------
if nargin < 3
    order = 2;
end

Nv = length(vert);
Nf = length(face);

edge  = sparse(Nv,Nv);
EuclD = ones(3,1);
for i = 1:Nf;
    
    % 1st order: compute Euclidean distances
    %----------------------------------------------------------------------
    if order
       Diff  = [vert(face(i,[1 2 3]),:) - vert(face(i,[2 3 1]),:)];
       EuclD = sqrt( sum(Diff.^2, 2) );
    end

    edge(face(i,1),face(i,2)) = EuclD(1);
    edge(face(i,2),face(i,3)) = EuclD(2);
    edge(face(i,3),face(i,1)) = EuclD(3);

    edge(face(i,2),face(i,1)) = EuclD(1);
    edge(face(i,3),face(i,2)) = EuclD(2);
    edge(face(i,1),face(i,3)) = EuclD(3);
end

% 0th order connectivity (Adjacency matrix)
%--------------------------------------------------------------------------
Mdist = edge;
if ~order
    Mdist = ~~Mdist;
end

% 2nd order connectivity
%--------------------------------------------------------------------------
if order == 2;
    for i = 1:Nv
        a          = find(edge(i,:));
        [b,c]      = find(edge(a,:));
        Mdist(i,c) = Mdist(i,a(b)) + diag(Mdist(a(b),c))';
        Mdist(c,i) = Mdist(i,c)';
    end
end

return
%==========================================================================