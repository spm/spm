function tests = test_spm_opm_hfc
% Unit Tests for spm_opm_hfc
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_opm_hfc_1(testCase)

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
S.li=1;
S.reg = 1;
H1 = spm_opm_vslm(S);

% create homogeneous fields 
beta =randn(3,size(D,2));
D(:,:,1)= H1*beta;
D.save();

% try and remove fields 
S=[];
S.D=D;
hD = spm_opm_hfc(S);

% check that shielding factor is greater than 100
Y = std(squeeze(D(:,:,1)),[],2);
res = std(squeeze(hD(:,:,1)),[],2);
ratio = mean(Y./res);

testCase.verifyTrue(ratio>100);
