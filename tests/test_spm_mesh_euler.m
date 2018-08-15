function tests = test_spm_mesh_euler
% Unit Tests for spm_mesh_euler
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_mesh_euler.m 7397 2018-08-15 11:04:26Z guillaume $

tests = functiontests(localfunctions);


function test_spm_mesh_euler_polyhedron(testCase)
M = spm_mesh_polyhedron('tetrahedron');
exp = 2;
act = spm_mesh_euler(M);
testCase.verifyEqual(act, exp);

M = spm_mesh_polyhedron('octahedron');
exp = 2;
act = spm_mesh_euler(M);
testCase.verifyEqual(act, exp);

M = spm_mesh_polyhedron('icosahedron');
exp = 2;
act = spm_mesh_euler(M);
testCase.verifyEqual(act, exp);

function test_spm_mesh_euler_sphere(testCase)

M = spm_mesh_sphere(4);
exp = 2;
act = spm_mesh_euler(M);
testCase.verifyEqual(act, exp);

M1 = spm_mesh_sphere(4);
M2 = spm_mesh_sphere(5);
M2.vertices = M2.vertices + 2;
M = spm_mesh_join([M1,M2]);
exp = 4;
act = spm_mesh_euler(M);
testCase.verifyEqual(act, exp);
