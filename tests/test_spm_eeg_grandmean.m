function tests = test_spm_eeg_grandmean
% Unit Tests for spm_eeg_merge
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_eeg_grandmean_1(testCase)

spm('defaults','eeg');

fname = fullfile(spm('Dir'),'tests','data','OPM','test_opm.mat');
tname = fullfile(spm('Dir'),'tests','data','OPM','test_opm2.dat');

D     = spm_eeg_load(fname);
D = chantype(D,1:110,'MEG');
D=type(D,'evoked');
D(:,:,:)=3;
D.save();

D2 = clone(D,tname);
D2(:,:,:)=6;
D2.save();


S=[];
S.D={D,D2};
Do = spm_eeg_grandmean(S);
mu = mean(Do(:));
D2.delete();
Do.delete();
D=type(D,'continuous');
D.save()
err = abs(mu-4.5);

testCase.verifyTrue(err < 1e-7);
