function tests = test_spm_mesh_refine
% Unit Tests for spm_mesh_refine
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_mesh_refine.m 7395 2018-08-14 14:09:21Z guillaume $

tests = functiontests(localfunctions);


function test_spm_mesh_refine_polyhedron(testCase)
P = spm_mesh_polyhedron('tetrahedron');
M = spm_mesh_refine(P);
exp = 4*size(P.faces,1);
act = size(M.faces,1);
testCase.verifyEqual(act, exp);

function test_spm_mesh_refine_gifti(testCase)
P = gifti(fullfile(spm('Dir'),'canonical','cortex_8196.surf.gii'));
P = export(P,'patch');
M = spm_mesh_refine(P);
exp = 4*size(P.faces,1);
act = size(M.faces,1);
testCase.verifyEqual(act, exp);

P.cdata = spm_mesh_curvature(P);
M = spm_mesh_refine(P);
testCase.verifyTrue(isfield(M,'cdata'));
exp = size(M.vertices,1);
act = size(M.cdata,1);
testCase.verifyEqual(act, exp);
exp = size(P.cdata,2);
act = size(M.cdata,2);
testCase.verifyEqual(act, exp);
