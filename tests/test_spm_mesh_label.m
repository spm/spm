function tests = test_spm_mesh_label
% Unit Tests for spm_mesh_label
%__________________________________________________________________________

% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_mesh_label_vertices(testCase)
P = export(gifti(fullfile(spm('Dir'),'canonical','cortex_8196.surf.gii')),'patch');
[C, N] = spm_mesh_label(P,'vertices');
testCase.verifyEqual(N, size(P.vertices,1)./[2 2]);
act = unique(C);
exp = [1;2];
testCase.verifyEqual(act, exp);
act = size(C);
exp = [sum(N) 1];
testCase.verifyEqual(act, exp);


function test_spm_mesh_label_faces(testCase)
P = export(gifti(fullfile(spm('Dir'),'canonical','cortex_8196.surf.gii')),'patch');
[C, N] = spm_mesh_label(P);
[C2, N2] = spm_mesh_label(P,'faces');
testCase.verifyEqual(C, C2);
testCase.verifyEqual(N, N2);
testCase.verifyEqual(N, size(P.vertices,1)./[2 2]);
act = unique(C);
exp = [1;2];
testCase.verifyEqual(act, exp);
act = size(C);
exp = [size(P.faces,1) 1];
testCase.verifyEqual(act, exp);
