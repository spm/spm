function N = spm_mesh_normals(M, unit)
% Compute (unit) normals of a surface mesh
% FORMAT N = spm_mesh_normals(M)
% M        - a patch structure or a handle to a patch
% unit     - boolean to indicate unit normals or not [default: false]
%
% N        - a [Nx3] array of (unit) normals on vertices
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_normals.m 5903 2014-03-03 14:56:10Z guillaume $

if nargin < 2, unit = false; end

if ~ishandle(M)
    f = figure('visible','off');
    M = patch(M, 'parent',axes('parent',f), 'visible', 'off');
end

N = double(get(M,'VertexNormals'));

if isempty(N)
    t = triangulation(double(get(M,'Faces')),double(get(M,'Vertices')));
    N = -double(t.vertexNormal);
end

try, close(f); end

if unit
    normN = sqrt(sum(N.^2,2));
    normN(normN < eps) = 1;
    N     = N ./ repmat(normN,1,3);
end