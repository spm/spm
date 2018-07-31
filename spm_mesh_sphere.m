function M = spm_mesh_sphere(N,M)
% Return a triangle mesh of a unit sphere
% N        - number of subdivision iterations [Default: 5]
% M        - initial triangle mesh [Default: 'icosahedron']
%
% M        - patch structure
%__________________________________________________________________________
%
% Computed using geodesic subdivisions of an icosahedron.
% See https://www.wikipedia.org/wiki/Geodesic_polyhedron
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id$


%-Check input arguments
%--------------------------------------------------------------------------
if nargin < 1, N = 5; end
if nargin < 2, M = 'icosahedron'; end
if ischar(M),  M = spm_mesh_polyhedron(M); end


%-Geodesic subdivisions
%==========================================================================
for i=1:N
    
    Nv = size(M.vertices,1);
    Nf = size(M.faces,1);
    
    F  = zeros(4*Nf,3);
    A  = spm_mesh_adjacency(M);
        
    for f=1:Nf
        
        T0 = M.faces(f,:);
        T1 = T0([2 3 1]);
        V  = (M.vertices(T0,:) + M.vertices(T1,:)) / 2;
        
        s = 1:3;
        b = [false false false];
        for j=1:3
            if A(T0(j),T1(j)) == 1
                s(j) = Nv + 1;
                A(T0(j),T1(j)) = s(j);
                A(T1(j),T0(j)) = s(j);
                Nv = s(j);
                b(j) = true;
            else
                s(j) = A(T0(j),T1(j));
            end
        end
        M.vertices(end+1:end+nnz(b),:) = V(b,:);
        T0(4:6) = s;
        
        F(4*f+(-3:0),:) = T0([1 4 6;4 2 5;6 5 3;4 5 6]);
        
    end

    M.faces = F;
    
end

%-Project vertices to unit sphere
%--------------------------------------------------------------------------
M.vertices = bsxfun(@rdivide,M.vertices,sqrt(sum(M.vertices.^2,2)));
