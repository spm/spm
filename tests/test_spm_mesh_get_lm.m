function tests = test_spm_mesh_get_lm
% Unit Tests for spm_mesh_get_lm
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_mesh_get_lm.m 6867 2016-09-12 15:04:44Z guillaume $

tests = functiontests(localfunctions);


function test_spm_mesh_get_lm_octahedron(testCase)
M = spm_mesh_polyhedron('octahedron');
T = [2 1 NaN 1 1 2]';

exp = [1 6];
act = spm_mesh_get_lm(M,T);
testCase.verifyEqual(act, exp);

function test_spm_mesh_get_lm_one(testCase)
M = spm_mesh_polyhedron('octahedron');
T = NaN(size(M.vertices,1),1);
T(3) = 1;

exp = 3;
act = spm_mesh_get_lm(M,T);
testCase.verifyEqual(act, exp);

function test_spm_mesh_get_lm_NaN(testCase)
M = spm_mesh_polyhedron('octahedron');
T = NaN(size(M.vertices,1),1);

exp = [];
act = spm_mesh_get_lm(M,T);
testCase.verifyEqual(act, exp);

function test_spm_mesh_get_lm_several(testCase)
M = spm_mesh_polyhedron('icosahedron');
T = [1 NaN NaN NaN NaN NaN 7 8 9 10 11 12]';

exp = [1 12];
act = spm_mesh_get_lm(M,T);
testCase.verifyEqual(act, exp);
