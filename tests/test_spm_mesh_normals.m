function tests = test_spm_mesh_normals
% Unit Tests for spm_mesh_normals
%__________________________________________________________________________

% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_mesh_normals_1(testCase)
M = spm_mesh_polyhedron('tetrahedron');

t = triangulation(M.faces,M.vertices);
exp = -double(t.vertexNormal);
act = spm_mesh_normals(M, true);
testCase.verifyEqual(act, exp ,'AbsTol', 10*eps);
