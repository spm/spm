function tests = test_spm_Ncdf
% Unit Tests for spm_Ncdf
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_Ncdf.m 7547 2019-03-18 17:19:38Z guillaume $

tests = functiontests(localfunctions);


function test_spm_Ncdf_1(testCase)
exp = 0.5; % normcdf(0,0,1)
act = spm_Ncdf(0,0,1);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_Ncdf_2(testCase)
exp = 0.5; % normcdf(0,0,4)
act = spm_Ncdf(0,0,4);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_Ncdf_3(testCase)
exp = 0.158655253931457; % normcdf(0,2,sqrt(4))
act = spm_Ncdf(0,2,4);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_Ncdf_4(testCase)
exp = 0.001349898031630; % normcdf(-4,2,sqrt(4))
act = spm_Ncdf(-4,2,4);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_Ncdf_5(testCase)
exp = 0.933192798731142; % normcdf(5,2,sqrt(4))
act = spm_Ncdf(5,2,4);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_Ncdf_6(testCase)
exp = 0.999968328758167; % normcdf(10,2,sqrt(4))
act = spm_Ncdf(10,2,4);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_Ncdf_7(testCase)
exp = [0.000173309675567  % normcdf(-4:10,4,sqrt(5))
   0.000872559349764
   0.003645179045768
   0.012673659338734
   0.036819135060151
   0.089856247439500
   0.185546684761349
   0.327360423009289
   0.500000000000000
   0.672639576990711
   0.814453315238651
   0.910143752560500
   0.963180864939849
   0.987326340661266
   0.996354820954232]';
act = spm_Ncdf(-4:10,4,5);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);
