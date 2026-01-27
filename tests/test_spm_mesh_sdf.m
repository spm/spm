classdef test_spm_mesh_sdf < matlab.unittest.TestCase
% Unit Tests for spm_mesh_sdf
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging



methods (Test)


function test_spm_mesh_sdf_(testCase)
M = gifti(fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii'));
M = export(M,'patch');
M.faces = double(M.faces);

V = spm_vol(fullfile(spm('Dir'),'canonical','single_subj_T1.nii'));

F = spm_mesh_sdf(M, V);

testCase.verifyClass(F,'double');
testCase.verifySize(F,V.dim);
testCase.verifyFalse(any(isnan(F(:))));
testCase.verifyTrue(any(F(:)>0));
testCase.verifyTrue(any(F(:)<0));
end
end % methods (Test)

end % classdef