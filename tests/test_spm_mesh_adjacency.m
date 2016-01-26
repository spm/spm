function tests = test_spm_mesh_adjacency
% Unit Tests for spm_mesh_adjacency
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_mesh_adjacency.m 6694 2016-01-26 17:09:11Z guillaume $

tests = functiontests(localfunctions);


function test_spm_mesh_adjacency_tetrahedron(testCase)
M = spm_mesh_polyhedron('tetrahedron');

exp = sparse(ones(4)-eye(4));
act = spm_mesh_adjacency(M);
testCase.verifyEqual(act, exp);

act = spm_mesh_adjacency(M.faces);
testCase.verifyEqual(act, exp);


function test_spm_mesh_adjacency_octahedron(testCase)
M = spm_mesh_polyhedron('octahedron');

exp = sparse([...
    0 1 1 1 1 0
    1 0 1 0 1 1
    1 1 0 1 0 1
    1 0 1 0 1 1
    1 1 0 1 0 1
    0 1 1 1 1 0 ]);
act = spm_mesh_adjacency(M);
testCase.verifyEqual(act, exp);

act = spm_mesh_adjacency(M.faces);
testCase.verifyEqual(act, exp);


function test_spm_mesh_adjacency_icosahedron(testCase)
M = spm_mesh_polyhedron('icosahedron');

exp = sparse([...
    0 1 1 1 1 1 0 0 0 0 0 0
    1 0 1 0 0 1 1 1 0 0 0 0
    1 1 0 1 0 0 0 1 1 0 0 0
    1 0 1 0 1 0 0 0 1 1 0 0
    1 0 0 1 0 1 0 0 0 1 1 0
    1 1 0 0 1 0 1 0 0 0 1 0
    0 1 0 0 0 1 0 1 0 0 1 1
    0 1 1 0 0 0 1 0 1 0 0 1
    0 0 1 1 0 0 0 1 0 1 0 1
    0 0 0 1 1 0 0 0 1 0 1 1
    0 0 0 0 1 1 1 0 0 1 0 1
    0 0 0 0 0 0 1 1 1 1 1 0 ]);
act = spm_mesh_adjacency(M);
testCase.verifyEqual(act, exp);

act = spm_mesh_adjacency(M.faces);
testCase.verifyEqual(act, exp);
