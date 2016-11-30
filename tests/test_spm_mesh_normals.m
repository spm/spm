function tests = test_spm_mesh_normals
% Unit Tests for spm_mesh_normals
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_mesh_normals.m 6955 2016-11-30 11:40:52Z guillaume $

tests = functiontests(localfunctions);


function test_spm_mesh_normals_1(testCase)
M = spm_mesh_polyhedron('tetrahedron');

t = triangulation(M.faces,M.vertices);
exp = -double(t.vertexNormal);
act = spm_mesh_normals(M, true);
testCase.verifyEqual(act, exp ,'AbsTol', 10*eps);
