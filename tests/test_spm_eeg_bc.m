function tests = test_spm_eeg_bc
% Unit Tests for spm_eeg_bc
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_eeg_bc_1(testCase)

spm('defaults','eeg');

fname = fullfile(spm('Dir'),'tests','data','OPM','test_opm.mat');
D     = spm_eeg_load(fname);
D = chantype(D,1:110,'MEG');
D.save();

D(:,1:500,1)= 10;
D(:,501:1000,1)= 20;
D.save();

S=[];
S.D=D;
S.timewin = [0 .5];
bD = spm_eeg_bc(S);

baseCorrected = mean(mean(bD(:,1:500,:)))==0;
actCorrected = mean(mean(bD(:,501:1000,:)))==10;


testCase.verifyTrue(baseCorrected & actCorrected);
