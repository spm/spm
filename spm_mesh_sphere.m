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

% Guillaume Flandin
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


%-Check input arguments
%--------------------------------------------------------------------------
if nargin < 1, N = 5; end
if nargin < 2, M = 'icosahedron'; end
if ischar(M),  M = spm_mesh_polyhedron(M); end


%-Geodesic subdivisions
%==========================================================================
for i=1:N
    
    M = spm_mesh_refine(M);
    
end

%-Project vertices to unit sphere
%--------------------------------------------------------------------------
M.vertices = bsxfun(@rdivide,M.vertices,sqrt(sum(M.vertices.^2,2)));
