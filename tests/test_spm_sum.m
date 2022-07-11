function tests = test_spm_sum
% Unit Tests for spm_sum
%__________________________________________________________________________

% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_sum_1(testCase)
X = magic(5);
exp = sum(X);
act = spm_sum(X);
testCase.verifyEqual(act, exp);
exp = sum(X,1);
act = spm_sum(X,1);
testCase.verifyEqual(act, exp);
exp = sum(X,2);
act = spm_sum(X,2);
testCase.verifyEqual(act, exp);
exp = 25*26/2;
act = spm_sum(X,[1 2]);
testCase.verifyEqual(act, exp);

function test_spm_sum_2(testCase)
X = 1:10;
exp = sum(X);
act = spm_sum(X);
testCase.verifyEqual(act, exp);
exp = sum(X,1);
act = spm_sum(X,1);
testCase.verifyEqual(act, exp);
exp = sum(X,2);
act = spm_sum(X,2);
testCase.verifyEqual(act, exp);
exp = sum(X(:));
act = spm_sum(X,[1 2]);
testCase.verifyEqual(act, exp);
exp = sum(X(:));
act = spm_sum(X,[1 2 3]);
testCase.verifyEqual(act, exp);
exp = sum(X(:));
act = spm_sum(X,[2 1]);
testCase.verifyEqual(act, exp);
exp = sum(X(:));
act = spm_sum(X,[3 1 2]);
testCase.verifyEqual(act, exp);

function test_spm_sum_3(testCase)
X = (1:10)';
exp = sum(X);
act = spm_sum(X);
testCase.verifyEqual(act, exp);
exp = sum(X,1);
act = spm_sum(X,1);
testCase.verifyEqual(act, exp);
exp = sum(X,2);
act = spm_sum(X,2);
testCase.verifyEqual(act, exp);
exp = sum(X(:));
act = spm_sum(X,[1 2]);
testCase.verifyEqual(act, exp);
exp = sum(X(:));
act = spm_sum(X,[1 2 3]);
testCase.verifyEqual(act, exp);
exp = sum(X(:));
act = spm_sum(X,[2 1]);
testCase.verifyEqual(act, exp);
exp = sum(X(:));
act = spm_sum(X,[3 1 2]);
testCase.verifyEqual(act, exp);

function test_spm_sum_4(testCase)
X = [magic(5);1:5];
X = cat(3,X,fliplr(X),flipud(X));
exp = sum(X);
act = spm_sum(X);
testCase.verifyEqual(act, exp);
exp = sum(X,1);
act = spm_sum(X,1);
testCase.verifyEqual(act, exp);
exp = sum(X,2);
act = spm_sum(X,2);
testCase.verifyEqual(act, exp);
exp = sum(X,3);
act = spm_sum(X,3);
testCase.verifyEqual(act, exp);
exp = sum(sum(X,1),2);
act = spm_sum(X,[1 2]);
testCase.verifyEqual(act, exp);
exp = sum(sum(X,2),3);
act = spm_sum(X,[2 3]);
testCase.verifyEqual(act, exp);
exp = sum(sum(X,3),1);
act = spm_sum(X,[3 1]);
testCase.verifyEqual(act, exp);
exp = sum(X(:));
act = spm_sum(X,[1 2 3]);
testCase.verifyEqual(act, exp);
exp = spm_sum(X,[1 2 3]);
act = spm_sum(X,[3 1 2]);
testCase.verifyEqual(act, exp);
exp = spm_sum(X,[1 2 3]);
act = spm_sum(X,[2 1 3]);
testCase.verifyEqual(act, exp);
