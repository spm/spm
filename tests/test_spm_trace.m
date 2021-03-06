function tests = test_spm_trace
% Unit Tests for spm_trace
%__________________________________________________________________________

% Copyright (C) 2015-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_trace_trace(testCase)
A = rand(3);
B = rand(3);

exp = trace(A*B);
act = spm_trace(A,B);
tol = 1e-6;
testCase.verifyEqual(act, exp,'AbsTol',tol);


function test_spm_trace_sum(testCase)
A = rand(3);
B = rand(3);

exp = sum(sum(A'.*B));
act = spm_trace(A,B);
tol = 1e-6;
testCase.verifyEqual(act, exp,'AbsTol',tol);


function test_spm_trace_nonsquare_1(testCase)
A = rand(3,5);
B = rand(5,3);

exp = trace(A*B);
act = spm_trace(A,B);
tol = 1e-6;
testCase.verifyEqual(act, exp,'AbsTol',tol);


function test_spm_trace_nonsquare_2(testCase)
A = rand(5,3);
B = rand(3,5);

exp = trace(A*B);
act = spm_trace(A,B);
tol = 1e-6;
testCase.verifyEqual(act, exp,'AbsTol',tol);
