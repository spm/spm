function tests = test_spm_invNcdf
% Unit Tests for spm_invNcdf
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_invNcdf.m 6860 2016-08-25 12:00:10Z guillaume $

tests = functiontests(localfunctions);


function test_spm_invNcdf_1(testCase)
exp = 1.644853626951473; % norminv(1 - 0.05)
act = spm_invNcdf(1 - 0.05);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_invNcdf_2(testCase)
exp = -1.644853626951473; % norminv(0.05)
act = spm_invNcdf(0.05);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_invNcdf_3(testCase)
exp = [NaN -Inf 0 Inf NaN];
ws = warning('off','SPM:outOfRangeNormal');
act = spm_invNcdf([-1 0 0.5 1 2]);
warning(ws);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_invNcdf_4(testCase)
exp = 4.848970052893895; % norminv(1 - 0.05,2,sqrt(3))
act = spm_invNcdf(1 - 0.05, 2, 3);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_invNcdf_5(testCase)
exp = -8.222082216130437; % norminv(10^-16)
act = spm_invNcdf(10^-16);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_invNcdf_6(testCase)
exp = -8.493793224109599; % norminv(10^-17)
act = spm_invNcdf(10^-17);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);
