function tests = test_spm_opm_create
% Unit Tests for spm_opm_hfc
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_opm_create_1(testCase)

spm('defaults','eeg');

data = fullfile(spm('Dir'),'tests','data','OPM','OPM_meg_001.cMEG');
pos = fullfile(spm('Dir'),'tests','data','OPM','OPM_HelmConfig.tsv');

S = [];
S.data = data;
S.positions = pos;
D = spm_opm_create(S);
test = -1.8334e4;
act = round(D(111,5));

testCase.verifyTrue((test-act)< 1e-6);


