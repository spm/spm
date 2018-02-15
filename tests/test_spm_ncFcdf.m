function tests = test_spm_ncFcdf
% Unit Tests for spm_ncFcdf
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_ncFcdf.m 7259 2018-02-15 10:16:32Z guillaume $

tests = functiontests(localfunctions);


function test_spm_ncFcdf_1(testCase)
exp = 0.659106867697940; % ncfcdf(1,1,10,0) = fcdf(1,1,10)
act = spm_ncFcdf(1,[1,10],0);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

exp = 0.603430543396035; % ncfcdf(1,2,12,0)
act = spm_ncFcdf(1,[2,12],0);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

exp = 0.571988250151978; % ncfcdf(2,2,10,2)
act = spm_ncFcdf(2,[2,10],2);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);

function test_spm_ncFcdf_2(testCase)
exp = [0 % ncfcdf(0:0.3:5,3,24,pi)
   0.048291661102670
   0.131225616302767
   0.226887358193830
   0.324414361264352
   0.417503323903781
   0.502767417050265
   0.578716229533475
   0.645046167955536
   0.702151693127709
   0.750799339445852
   0.791918495436190
   0.826473124442134
   0.855387792826706
   0.879508986629912
   0.899588574806357
   0.916280590482256]';
act = spm_ncFcdf(0:0.3:5,[3,24],pi);
tol = 1e-12;
testCase.verifyEqual(act, exp,'AbsTol',tol);
