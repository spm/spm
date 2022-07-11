function tests = test_spm_gamrnd
% Unit Tests for spm_gamrnd
%__________________________________________________________________________

% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_gamrnd_1(testCase)
import matlab.unittest.constraints.HasSize
import matlab.unittest.constraints.IsGreaterThan

g = spm_gamrnd(2,5);
testCase.verifyThat(g, HasSize([1 1]));
testCase.verifyThat(g, IsGreaterThan(0));

g = spm_gamrnd(2,5,2);
%testCase.verifyThat(g, HasSize([2 2]));
testCase.verifyThat(g, IsGreaterThan(0));

g = spm_gamrnd(2,5,3,5);
testCase.verifyThat(g, HasSize([3 5]));
testCase.verifyThat(g, IsGreaterThan(0));

g = spm_gamrnd(2,5,[3,5]);
%testCase.verifyThat(g, HasSize([3 5]));
testCase.verifyThat(g, IsGreaterThan(0));
