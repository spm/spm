function tests = test_spm_get_data
% Unit Tests for spm_get_data
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_get_data.m 7385 2018-08-01 15:58:24Z guillaume $

tests = functiontests(localfunctions);


function test_spm_get_data_1(testCase)
act = spm_get_data(spm_vol(fullfile(spm('Dir'),'tpm','mask_ICV.nii')),[1 1 1]');
testCase.verifyTrue(isnumeric(act));
testCase.verifyTrue(isequal(numel(act),1));

act = spm_get_data(fullfile(spm('Dir'),'tpm','mask_ICV.nii'),[1 1 1;2 2 2]');
testCase.verifyTrue(isnumeric(act));
testCase.verifyTrue(isequal(size(act),[1 2]));

act = spm_get_data(char(...
    fullfile(spm('Dir'),'canonical','avg152PD.nii'),...
    fullfile(spm('Dir'),'canonical','avg152T1.nii'),...
    fullfile(spm('Dir'),'canonical','avg152T2.nii')),[1 1 1;2 2 2]');
testCase.verifyTrue(isnumeric(act));
testCase.verifyTrue(isequal(size(act),[3 2]));


function test_spm_get_data_error(testCase)
err = @() evalc('spm_get_data(struct(''fname'',tempname),[1 1 1]'')');
testCase.verifyError(err, ?MException);
