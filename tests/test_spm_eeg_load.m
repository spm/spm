function tests = test_spm_eeg_load
% Unit Tests for spm_eeg_load
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_eeg_load_1(testCase)

spm('defaults','eeg');

fname = fullfile(spm('Dir'),'tests','data','OPM','test_opm.mat');
D     = spm_eeg_load(fname);

exp = 'meeg';
act = class(D);
testCase.verifyTrue(isequal(exp, act));

exp = 1;
act = numel(D);
testCase.verifyTrue(isequal(exp, act));

exp = [110,1000,1];
act = size(D);
testCase.verifyTrue(isequal(exp, act));
