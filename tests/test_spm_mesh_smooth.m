function tests = test_spm_mesh_smooth
% Unit Tests for spm_mesh_smooth
%__________________________________________________________________________

% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_mesh_smooth_1(testCase)

M = spm_mesh_polyhedron('tetrahedron');
K = spm_mesh_smooth(M);
exp = [size(M.vertices,1) size(M.vertices,1)];
act = size(K);
testCase.verifyEqual(act, exp);

exp = K;
act = spm_mesh_smooth(K);
testCase.verifyEqual(act, exp);


function test_spm_mesh_smooth_2(testCase)

M = spm_mesh_sphere(4);
K = spm_mesh_smooth(M);
T = rand(size(M.vertices,1),1);
T = spm_mesh_smooth(K,T,10);
testCase.verifyTrue(min(T)>=0);
testCase.verifyTrue(max(T)<=1);
