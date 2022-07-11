function tests = test_spm_mesh_neighbours
% Unit Tests for spm_mesh_neighbours
%__________________________________________________________________________

% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_mesh_neighbours_tetrahedron(testCase)
M = spm_mesh_polyhedron('tetrahedron');

exp = [2 3 4;1 3 4;1 2 4;1 2 3];
act = spm_mesh_neighbours(M);
testCase.verifyEqual(act, exp);

act = spm_mesh_neighbours(sparse(ones(4)-eye(4)));
testCase.verifyEqual(act, exp);
