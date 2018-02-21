function tests = test_spm_ncFpdf
% Unit Tests for spm_ncFpdf
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_ncFpdf.m 7263 2018-02-21 13:30:16Z guillaume $

tests = functiontests(localfunctions);


function test_spm_ncFpdf_1(testCase)
exp = 0.230361989229139; % ncfpdf(1,1,10,0) = fpdf(1,1,10)
act = spm_ncFpdf(1,[1,10],0);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

exp = 0.339916677089114; % ncfpdf(1,2,12,0)
act = spm_ncFpdf(1,[2,12],0);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

exp = 0.187052298433871; % ncfpdf(2,2,10,2)
act = spm_ncFpdf(2,[2,10],2);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_ncFpdf_2(testCase)
exp = [0 % ncfpdf(0:0.3:5,3,24,pi)
   0.237214867508572
   0.305378558226511
   0.326510412867316
   0.320250126742450
   0.298516612369365
   0.269121488384202
   0.237061968288903
   0.205374511459040
   0.175766182458195
   0.149074704100068
   0.125592498668479
   0.105287389552812
   0.087948186675270
   0.073277611204176
   0.060949393292047
   0.050641638818210]';
act = spm_ncFpdf(0:0.3:5,[3,24],pi);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);
