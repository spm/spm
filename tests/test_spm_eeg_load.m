function tests = test_spm_eeg_load
% Unit Tests for spm_eeg_load
%__________________________________________________________________________
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging

tests = functiontests(localfunctions);


function test_spm_eeg_load_1(testCase)
fname = fullfile(spm('dir'),'tests','data','test_opm.mat');
D = spm_eeg_load(fname);
act = size(D);
exp=[110,1000,1];
testCase.verifyEqual(act, exp);

