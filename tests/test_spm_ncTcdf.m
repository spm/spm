function tests = test_spm_ncTcdf
% Unit Tests for spm_ncTcdf
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_ncTcdf.m 7258 2018-02-14 13:09:46Z guillaume $

tests = functiontests(localfunctions);


function test_spm_ncTcdf_1(testCase)
exp = 0.5; % nctcdf(0,1,0)
act = spm_ncTcdf(0,1,0);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_ncTcdf_2(testCase)
exp = 0.5; % nctcdf(0,4,0)
act = spm_ncTcdf(0,4,0);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_ncTcdf_3(testCase)
exp = 0.022750131948179; % nctcdf(0,4,2)
act = spm_ncTcdf(0,4,2);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_ncTcdf_4(testCase)
exp = 4.218150123125319e-05; % nctcdf(-4,4,2)
act = spm_ncTcdf(-4,4,2);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_ncTcdf_5(testCase)
exp = 0.920822617476195; % nctcdf(5,4,2)
act = spm_ncTcdf(5,4,2);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_ncTcdf_6(testCase)
exp = 0.992600181884681; % nctcdf(10,4,2)
act = spm_ncTcdf(10,4,2);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_ncTcdf_7(testCase)
exp = [0.000000011096548  % nctcdf(-4:10,4,4)
   0.000000032879594
   0.000000140813044
   0.000001174173412
   0.000031671241833
   0.001992604459483
   0.042091256308006
   0.202851164067265
   0.429784541497594
   0.622224194923879
   0.754499466251885
   0.838999243587669
   0.892202382905940
   0.926042984950848
   0.948002154097108]';
act = spm_ncTcdf(-4:10,4,4);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);
