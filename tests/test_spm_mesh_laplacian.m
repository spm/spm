function tests = test_spm_mesh_laplacian
% Unit Tests for spm_mesh_laplacian
%__________________________________________________________________________

% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_mesh_laplacian_sphere(testCase)
M = spm_mesh_sphere;

L  = spm_mesh_laplacian(M);
Lg = spm_mesh_laplacian(M,'graph');
Lm = spm_mesh_laplacian(M,'mesh');

testCase.verifyEqual(L, Lg);

exp = [size(M.vertices,1) size(M.vertices,1)];
act = size(Lg);
testCase.verifyEqual(act, exp);
act = size(Lm);
testCase.verifyEqual(act, exp);

exp = true;
act = issparse(Lg);
testCase.verifyEqual(act, exp);
act = issparse(Lm);
testCase.verifyEqual(act, exp);

testCase.verifyEqual(nnz(Lg), nnz(Lm));
testCase.verifyEqual(find(Lg), find(Lm));

exp = zeros(size(M.vertices,1),1);
act = full(sum(Lg,2));
testCase.verifyEqual(act, exp);
act = full(sum(Lm,2));
tol = 1e-10;
testCase.verifyEqual(act, exp,'AbsTol',tol);
