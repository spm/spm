function tests = test_spm_eeg_average
% Unit Tests for spm_eeg_average
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_eeg_average_1(testCase)

spm('defaults','eeg');

fname = fullfile(spm('Dir'),'tests','data','OPM','test_opm.mat');
D     = spm_eeg_load(fname);
D = chantype(D,1,'TRIG');
D.save();
D(:,:,:)=0;
D(1,30:60,1)=3;
D(1,200:230,1)=3;
D(1,700:730,1)=3;
D.save();


S=[];
S.D=D;
S.timewin= [-10,20];
eD=spm_opm_epoch_trigger(S);

eD(2:end,:,1)=3;
eD(2:end,:,2)=6;
eD(2:end,:,3)=9;

S=[];
S.D=eD;
mD = spm_eeg_average(S);

act = mean(mean(squeeze(mD(2:end,:,:))));
exp = 6;


testCase.verifyTrue(abs((act-exp)/exp)< 1e-6);
