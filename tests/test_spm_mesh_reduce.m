function tests = test_spm_mesh_reduce
% Unit Tests for spm_mesh_reduce
%__________________________________________________________________________

% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_mesh_sphere(testCase)
M = spm_mesh_sphere;
m = spm_mesh_reduce(M,1000);
testCase.verifyTrue(size(m.faces,1) <= 1000);

function test_spm_mesh_cortex(testCase)
M = gifti(fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii'));
M = export(M,'patch');
for i=round(linspace(1000,size(M.faces,1),50))
    m = spm_mesh_reduce(M,i);
    testCase.verifyTrue(size(m.faces,1) <= i);
end
