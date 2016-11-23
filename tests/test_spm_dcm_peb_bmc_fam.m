function tests = test_spm_dcm_peb_bmc_fam
% Unit Tests for test_spm_dcm_peb_bmc
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dcm_peb_bmc_fam.m 6946 2016-11-23 15:26:29Z peter $

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_peb_bmc(testCase)

data_path = get_data_path();

% Load
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
PEB = load(fullfile(data_path,'PEB_test.mat'));
GCM = GCM.GCM;
PEB = PEB.PEB;

% Model comparison
[BMA,BMR] = spm_dcm_peb_bmc(PEB(1),GCM(1,:));

% Family 1: model 1 (full) and model 2 (fwd modulation)
% Family 2: model 3 (bwd) and model 4 (bwd modulation)
families = [1 1 2 2];

% Model comparison under families
[BMA,fam] = spm_dcm_peb_bmc_fam(BMA,BMR,families);

% -------------------------------------------------------------------------
% Should be equal priors as families were balanced
Nm = size(GCM,2);
testCase.assertEqual([0.5 0.5], fam.family.prior)
testCase.assertTrue(all(fam.model.prior(:) == 1/(Nm*Nm)));

% Family 1 and model 2 should win
testCase.assertTrue(fam.family.post(1,1) > 0.9);
testCase.assertTrue(fam.model.Pw(2) > 0.9);
testCase.assertTrue(fam.model.Px(2) > 0.9);

% -------------------------------------------------------------------------
% Run again with unbalanced families
families = [1 1 1 2];
[BMA,fam] = spm_dcm_peb_bmc_fam(BMA,BMR,families);

% Families should have equal priors
testCase.assertEqual([0.5 0.5], fam.family.prior);

% Model joint prior matrix should sum to 1
testCase.assertEqual(1, sum(fam.model.prior(:)));

% Models within each family should have same total prior probability
family_total_prior = [];
for f = 1:2
    p = (fam.model.prior(families == f,families == f));
    family_total_prior(f) = sum(p(:));
end
testCase.assertTrue(...
    all(family_total_prior(:) - family_total_prior(1) < 1e-10));

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');