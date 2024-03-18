function tests = test_spm_opm_epoch_trigger
% Unit Tests for spm_eeg_average
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_opm_epoch_trigger_1(testCase)

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

exp = 3;
act = size(eD,3);
testCase.verifyTrue(isequal(exp, act));
