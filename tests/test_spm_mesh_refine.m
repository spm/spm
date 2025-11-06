classdef test_spm_mesh_refine < matlab.unittest.TestCase
% Unit Tests for spm_mesh_refine
%__________________________________________________________________________

% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging



methods (Test)


function test_spm_mesh_refine_polyhedron(testCase)
P = spm_mesh_polyhedron('tetrahedron');
M = spm_mesh_refine(P);
exp = 4*size(P.faces,1);
act = size(M.faces,1);
testCase.verifyEqual(act, exp);
end

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
end
end % methods (Test)

end % classdef