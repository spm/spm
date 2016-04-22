function tests = test_spm_z2p
% Unit Tests for spm_z2p
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_z2p.m 6778 2016-04-22 11:51:29Z guillaume $

tests = functiontests(localfunctions);


function test_spm_z2p_Z(testCase)
exp = 0.05;
act = spm_z2p(1.644853626951472,[],'Z');
tol = 1e-8;
testCase.verifyEqual(act, exp,'AbsTol',tol);


function test_spm_z2p_T(testCase)
exp = 0.05;
act = spm_z2p(1.745883689098006,[1 16],'T');
tol = 1e-8;
testCase.verifyEqual(act, exp,'AbsTol',tol);


function test_spm_z2p_X(testCase)
exp = 0.05;
act = spm_z2p(26.296227607876062,[1 16],'X');
tol = 1e-8;
testCase.verifyEqual(act, exp,'AbsTol',tol);


function test_spm_z2p_F(testCase)
exp = 0.05;
act = spm_z2p(3.633723519202212,[2 16],'F');
tol = 1e-8;
testCase.verifyEqual(act, exp,'AbsTol',tol);


function test_spm_z2p_P(testCase)
exp = 0.5;
act = spm_z2p(0.5,[],'P');
tol = 1e-8;
testCase.verifyEqual(act, exp,'AbsTol',tol);
