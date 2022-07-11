function tests = test_spm_mesh_ray_intersect
% Unit Tests for spm_mesh_ray_intersect
%__________________________________________________________________________

% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_mesh_ray(testCase)
M = spm_mesh_sphere(3);

R = struct('orig',[0.1 2 0.2]','vec',[0 -1 0]');
[I,P] = spm_mesh_ray_intersect(M,R);
testCase.verifyTrue(islogical(I));
exp = 2;
act = nnz(I);
testCase.verifyEqual(act, exp);
exp = size(M.faces,1);
act = numel(I);
testCase.verifyEqual(act, exp);
exp = [2 3];
act = size(P);
testCase.verifyEqual(act, exp);

R = struct('orig',[0.1 2 0.2]','vec',[1 0 0]');
[I,P] = spm_mesh_ray_intersect(M,R);
exp = 0;
act = nnz(I);
testCase.verifyEqual(act, exp);
exp = [0 3];
act = size(P);
testCase.verifyEqual(act, exp);
