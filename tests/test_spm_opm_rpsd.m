function tests = test_spm_opm_rpsd
% Unit Tests for spm_opm_rpsd
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_opm_rpsd_1(testCase)

spm('defaults','eeg');

fname = fullfile(spm('Dir'),'tests','data','OPM','test_opm.mat');
tname = fullfile(spm('Dir'),'tests','data','OPM','test_opm2.dat');

D     = spm_eeg_load(fname);
D = chantype(D,1:110,'MEG');
D.save();
D(:,:,:)=randn(size(D));
D.save();

D2 = clone(D,tname);
D2(:,:,:)=randn(size(D2))*3;
D2.save()

S=[];
S.D1=D;
S.D2=D2;
S.triallength = 1000;
S.dB = 1;
[shield,~] = spm_opm_rpsd(S);
act = mean(mean(shield));
D2.delete();

% test  for at least 9 dB of shielding
testCase.verifyTrue(act < -9);
