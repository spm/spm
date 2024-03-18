function tests = test_spm_opm_headmodel
% Unit Tests for spm_opm_headmodel
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_opm_headmodel_1(testCase)

spm('defaults','eeg');

fname = fullfile(spm('Dir'),'tests','data','OPM','test_opm.mat');
D     = spm_eeg_load(fname);
D = chantype(D,1:110,'MEG');
D.inv{1}.mesh.sMRI = fullfile(spm('Dir'),'canonical','single_subj_T1.nii');
D.inv{1}.mesh.tess_ctx = fullfile(spm('Dir'),'canonical','cortex_5124.surf.gii');
D.inv{1}.mesh.tess_scalp = fullfile(spm('Dir'),'canonical','scalp_2562.surf.gii');
D.inv{1}.mesh.tess_oskull = fullfile(spm('Dir'),'canonical','oskull_2562.surf.gii');
D.inv{1}.mesh.tess_iskull = fullfile(spm('Dir'),'canonical','iskull_2562.surf.gii');
D.save();

S=[];
S.D =D;
S.template = 1;
S.meshres= 1;
S.lead = 1;
[~,L] = spm_opm_headmodel(S);

C = full(L)*full(L)';
observed = sum((cumsum(svd(C))/sum(svd(C)))<.99);
expected = 50;

% check dimensionality > 50
testCase.verifyTrue(observed>=expected);
