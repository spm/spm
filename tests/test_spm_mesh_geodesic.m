function tests = test_spm_mesh_geodesic
% Unit Tests for spm_mesh_geodesic
%__________________________________________________________________________

% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_mesh_geodesic_(testCase)
M = spm_mesh_sphere;
D = spm_mesh_geodesic(M,1);
exp = pi;
act = max(D);
testCase.verifyEqual(act, exp, 'AbsTol',1e-3);

[~,L] = spm_mesh_geodesic(M,1);
exp = 1;
act = unique(L);
testCase.verifyEqual(act, exp);

[~,~,P] = spm_mesh_geodesic(M,1);
testCase.verifyTrue(iscell(P));

[D,L] = spm_mesh_geodesic(M,1:100:size(M.vertices,1));
testCase.verifyTrue(max(D) < pi);
testCase.verifyTrue(numel(unique(L)) > 1);
