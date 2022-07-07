function L = spm_mesh_laplacian(M,T)
% Compute the graph or (cotangent) mesh Laplacian
% M        - patch structure: vertices and faces must be mx3 and nx3 arrays
% T        - {'graph','mesh'} [Default: 'graph']
%
% L        - Laplacian
%__________________________________________________________________________
%
% Laplacian matrix:
%   https://en.wikipedia.org/wiki/Laplacian_matrix
%   https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Mesh_Laplacians
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


if nargin < 2, T = 'graph'; end

if ismember(T,{'graph','combinatorial','topological'})
    
    %-Adjacency matrix
    A = spm_mesh_distmtx(M,0);
    
    %-Degree matrix
    n = size(A,1);
    D = spdiags(sum(A,2),0,n,n);
    
    %-Graph Laplacian
    L = D - A;
    
elseif ismember(T,{'mesh','geometric','cotangent'})
    
    F = M.faces;
    
    %-Edge lengths
    l = M.vertices(F',:);
    l = permute(reshape(l',3,3,[]),[2 1 3]);
    l = squeeze(sqrt(sum((l([1 2 3],:,:) - l([2 3 1],:,:)).^2,2)));
    
    %-Face areas
    Af  = spm_mesh_area(l,'face');
    
    %-Cotan weights
    %cot = (1-2*eye(3)) * l.^2 ./ repmat(4*Af,3,1);
    cot = bsxfun(@rdivide, (1-2*eye(3)) * l.^2, 4*Af);
    
    %-Cotangent matrix
    i = [F(:,1);F(:,2);F(:,2);F(:,3);F(:,1);F(:,3);F(:,1);F(:,2);F(:,3)];
    j = [F(:,2);F(:,1);F(:,3);F(:,2);F(:,3);F(:,1);F(:,1);F(:,2);F(:,3)];
    s = [cot(1,:),cot(1,:),cot(2,:),cot(2,:),cot(3,:),cot(3,:),...
        -cot(1,:)-cot(3,:),-cot(1,:)-cot(2,:),-cot(2,:)-cot(3,:)] / 2;
    C = sparse(i,j,s,size(M.vertices,1),size(M.vertices,1));
    
    %-Mesh Laplacian (without mass matrix weighting)
    L = C;

else
    
    error('Unknown option.');
    
end
