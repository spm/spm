function tests = test_spm_ncTpdf
% Unit Tests for spm_ncTpdf
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_ncTpdf.m 7260 2018-02-19 10:55:53Z guillaume $

tests = functiontests(localfunctions);


function test_spm_ncTpdf_1(testCase)
exp = spm_Tpdf(0,1);
act = spm_ncTpdf(0,1,0);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

exp = spm_Tpdf(-2:0.5:2,2);
act = spm_ncTpdf(-2:0.5:2,2,0);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

exp = spm_Tpdf(0,1:4);
act = spm_ncTpdf(0,1:4,0);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

exp = spm_Tpdf(-2:2,1:5);
act = spm_ncTpdf(-2:2,1:5,0);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_ncTpdf_2(testCase)
exp = 0.193064705260108; % nctpdf(0,1,1)
act = spm_ncTpdf(0,1,1);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

exp = 2.072225240523640e-04; % nctpdf(1,10,5);
act = spm_ncTpdf(1,10,5);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

exp = [... % nctpdf(-1:3,6:10,-2:2);
   0.240403340723731
   0.233509118657019
   0.227607580145303
   0.223756811547761
   0.219731001980809]';
act = spm_ncTpdf(-1:3,6:10,-2:2);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);
