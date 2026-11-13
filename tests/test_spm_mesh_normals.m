classdef test_spm_mesh_normals < matlab.unittest.TestCase
% Unit Tests for spm_mesh_normals
%__________________________________________________________________________

% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging



methods (Test)


function test_spm_mesh_normals_1(testCase)
M = spm_mesh_polyhedron('tetrahedron');

t = triangulation(M.faces,M.vertices);
exp = -double(t.vertexNormal);
act = spm_mesh_normals(M, true);
testCase.verifyEqual(act, exp ,'AbsTol', 10*eps);
end
end % methods (Test)

end % classdef