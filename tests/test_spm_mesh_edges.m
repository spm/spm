classdef test_spm_mesh_edges < matlab.unittest.TestCase
% Unit Tests for spm_mesh_edges
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging



methods (Test)


function test_spm_mesh_edges_struct(testCase)
M = spm_mesh_cube;

E = spm_mesh_edges(M);

exp = [18, 2];
act = size(E);
testCase.verifyEqual(act, exp);

[E,L] = spm_mesh_edges(M);

exp = [18, 1];
act = size(L);
testCase.verifyEqual(act, exp);

%exp = 2;
%act = numel(uniquetol(L));
%testCase.verifyEqual(act, exp);

exp = true;
act = all(L < 2);
testCase.verifyEqual(act, exp);
end

function test_spm_mesh_edges_faces(testCase)
M = spm_mesh_cube;

E = spm_mesh_edges(M.faces);

exp = [18, 2];
act = size(E);
testCase.verifyEqual(act, exp);
end
end % methods (Test)

end % classdef