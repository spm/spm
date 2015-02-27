function tests = test_spm_dctmtx
% Unit Tests for spm_dctmtx
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dctmtx.m 6352 2015-02-27 18:30:35Z guillaume $

tests = functiontests(localfunctions);


function test_spm_dctmtx_1(testCase)
N = 16;
C = spm_dctmtx(N);

import matlab.unittest.constraints.HasSize
testCase.verifyThat(C, HasSize([N N]));

exp = eye(N);
act = C'*C;
tol = 100*eps;
testCase.verifyEqual(act, exp,'AbsTol',tol);


function test_spm_dctmtx_2(testCase)
N = 16;
K = 8;
C = spm_dctmtx(N,K);

exp = [N K];
act = size(C);
testCase.verifyEqual(act, exp);

exp = eye(K);
act = C'*C;
tol = 100*eps;
testCase.verifyEqual(act, exp,'AbsTol',tol);


function test_spm_dctmtx_3(testCase)
N = 16;
K = 8;
C = spm_dctmtx(N,K,1:2:N);

exp = [N/2 K];
act = size(C);
testCase.verifyEqual(act, exp);

C = spm_dctmtx(N,K,'diff');

exp = [N K];
act = size(C);
testCase.verifyEqual(act, exp);

C = spm_dctmtx(N,K,'diff2');

exp = [N K];
act = size(C);
testCase.verifyEqual(act, exp);


function test_spm_dctmtx_4(testCase)
N = 16;
K = 8;
n = 1:2:N;

C = spm_dctmtx(N,K,n,'diff');

exp = [N/2 K];
act = size(C);
testCase.verifyEqual(act, exp);

C = spm_dctmtx(N,K,n,'diff2');

exp = [N/2 K];
act = size(C);
testCase.verifyEqual(act, exp);

testCase.verifyError(@()spm_dctmtx(N,K,n,'XXX'),?MException);
