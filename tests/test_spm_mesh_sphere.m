function tests = test_spm_mesh_sphere
% Unit Tests for spm_mesh_sphere
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_mesh_sphere.m 7395 2018-08-14 14:09:21Z guillaume $

tests = functiontests(localfunctions);


function test_spm_mesh_sphere_icosahedron(testCase)
M = spm_mesh_sphere(0);
testCase.verifyTrue(isstruct(M));
exp = 20;
act = size(M.faces,1);
testCase.verifyEqual(act, exp);
exp = 12;
act = size(M.vertices,1);
testCase.verifyEqual(act, exp);

M = spm_mesh_sphere(1);
testCase.verifyTrue(isstruct(M));
exp = [42 80];
act = [size(M.vertices,1) size(M.faces,1)];
testCase.verifyEqual(act, exp);

M = spm_mesh_sphere(2);
testCase.verifyTrue(isstruct(M));
exp = [162 320];
act = [size(M.vertices,1) size(M.faces,1)];
testCase.verifyEqual(act, exp);

M = spm_mesh_sphere(3);
testCase.verifyTrue(isstruct(M));
exp = [642 1280];
act = [size(M.vertices,1) size(M.faces,1)];
testCase.verifyEqual(act, exp);

M = spm_mesh_sphere(4);
testCase.verifyTrue(isstruct(M));
exp = [2562 5120];
act = [size(M.vertices,1) size(M.faces,1)];
testCase.verifyEqual(act, exp);

M = spm_mesh_sphere(5,'icosahedron');
testCase.verifyTrue(isstruct(M));
exp = [10242 20480];
act = [size(M.vertices,1) size(M.faces,1)];
testCase.verifyEqual(act, exp);


function test_spm_mesh_sphere_octahedron(testCase)
M = spm_mesh_sphere(0,spm_mesh_polyhedron('octahedron'));
testCase.verifyTrue(isstruct(M));
exp = [6 8];
act = [size(M.vertices,1) size(M.faces,1)];
testCase.verifyEqual(act, exp);

M = spm_mesh_sphere(5,'octahedron');
testCase.verifyTrue(isstruct(M));
exp = [4098 8192];
act = [size(M.vertices,1) size(M.faces,1)];
testCase.verifyEqual(act, exp);
