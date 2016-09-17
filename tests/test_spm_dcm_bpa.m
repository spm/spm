function tests = test_spm_dcm_bpa
% Unit Tests for spm_cfg_dcm_peb (PEB batch)
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dcm_bpa.m 6880 2016-09-17 18:00:38Z peter $

tests = functiontests(localfunctions);

function test_bpa(testCase)
% A simple test in the absence of covariance

% Prepare DCM 1
DCM1      = struct();
DCM1.M.pC = 1;
DCM1.M.pE = 0;
DCM1.Cp   = 0.5;
DCM1.Ep   = -0.5;

% Prepare DCM 2
DCM2    = DCM1;
DCM2.Ep = 0.5;

% Collate
GCM = {DCM1;DCM2};

% Average
BPA = spm_dcm_bpa(GCM);
BPA.Cp = full(BPA.Cp);

% Check. Precisions add, means are precision-weighted and added
expected = 1 / (1/DCM1.Cp + 1/DCM2.Cp);
testCase.assertEqual(BPA.Ep, 0,'AbsTol',0.00001);
testCase.assertEqual(BPA.Cp, expected,'AbsTol',0.00001);

% Posterior probability should be 0.5
testCase.assertEqual(BPA.Pp, 0.5,'AbsTol',0.00001);

% -------------------------------------------------------------------------
function test_bpa_pruned(testCase)
% Test with a parameter fixed at zero

% Prepare DCM 1
DCM1      = struct();
DCM1.M.pC = diag([1 0]);
DCM1.M.pE = [0 0];
DCM1.Cp   = diag([0.5 0]);
DCM1.Ep   = [-0.5 0];

% Prepare DCM 2
DCM2    = DCM1;
DCM2.Ep = [0.5 0];

% Collate
GCM = {DCM1;DCM2};

% Average
BPA = spm_dcm_bpa(GCM);
Cp  = full(diag(BPA.Cp));

% Check. Precisions add, means are precision-weighted and added
expected = 1 / (1 / DCM1.Cp(1) + 1 / DCM2.Cp(1));
testCase.assertEqual(BPA.Ep, [0 0],'AbsTol',0.00001);
testCase.assertEqual(Cp, [expected 0]','AbsTol',0.00001);

testCase.assertEqual(BPA.Pp, [0.5 NaN],'AbsTol',0.00001);

% -------------------------------------------------------------------------
function test_bpa_pruned2(testCase)
% Test with a fixed (disabled) parameter which has non-zero prior mean

% Prepare DCM 1
DCM1      = struct();
DCM1.M.pC = diag([1 0]);
DCM1.M.pE = [0 0.1];
DCM1.Cp   = diag([0.2 0]);
DCM1.Ep   = [1 0.1];

% Prepare DCM 2
DCM2    = DCM1;
DCM2.Ep = [2 0.1];

% Collate
GCM = {DCM1;DCM2};

% Average
BPA = spm_dcm_bpa(GCM);
Cp  = full(diag(BPA.Cp));

% Check. Precisions add, means are precision-weighted and added
expected1 = 1 / (1 / DCM1.Cp(1) + 1 / DCM2.Cp(1));
expected2 = 1 / (1 / DCM1.Cp(2) + 1 / DCM2.Cp(2));
testCase.assertEqual(BPA.Ep, [1.5 0.1],'AbsTol',0.00001);
testCase.assertEqual(Cp, [expected1 expected2]','AbsTol',0.00001);

testCase.assertEqual(BPA.Pp, [1 NaN],'AbsTol',0.00001);