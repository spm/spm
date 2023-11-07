function tests = test_spm_opm_psd
% Unit Tests for spm_eeg_average
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_opm_psd_1(testCase)

spm('defaults','eeg');

fname = fullfile(spm('Dir'),'tests','data','OPM','test_opm.mat');
D     = spm_eeg_load(fname);
D = chantype(D,1:110,'MEG');
D.save();
D(:,:,:)=randn(size(D));
D.save();


S=[];
S.D=D;
S.triallength=800;
S.plot=0;
[p,~]=spm_opm_psd(S);
act = mean(mean(p))*sqrt(800);

testCase.verifyTrue(act>.9);
