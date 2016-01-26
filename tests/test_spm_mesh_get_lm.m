function tests = test_spm_mesh_get_lm
% Unit Tests for spm_mesh_get_lm
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_mesh_get_lm.m 6694 2016-01-26 17:09:11Z guillaume $

tests = functiontests(localfunctions);


function test_spm_mesh_get_lm_octahedron(testCase)
M = spm_mesh_polyhedron('octahedron');
T = [2 1 NaN 1 1 2]';

exp = [1 6];
act = spm_mesh_get_lm(M,T);
testCase.verifyEqual(act, exp);
