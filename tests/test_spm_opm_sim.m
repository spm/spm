function tests = test_spm_opm_sim
% Unit Tests for spm_opm_sim
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_opm_sim_1(testCase)

spm('defaults','eeg');


S= [];
S.wholehead = 1;
S.nDens = 1;
S.lead = 1;
[D,L] = spm_opm_sim(S);
close all

exp = 'meeg';
act = class(D);
testCase.verifyTrue(isequal(exp, act));

exp = 1;
act = numel(D);
testCase.verifyTrue(isequal(exp, act));

exp = size(D,1);
act = size(L,1);
testCase.verifyTrue(isequal(exp, act));




