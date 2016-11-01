function tests = test_spm_dcm_bmr_all
% Unit Tests for test_spm_dcm_bmr_all
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dcm_bmr_all.m 6915 2016-11-01 17:38:08Z peter $

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_bmr_all(testCase)

data_path = get_data_path();

% Load PEB of full DCM
PEB = load(fullfile(data_path,'PEB_test.mat'));
PEB = PEB.PEB;
PEB = PEB(1);

% Load DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;

% Prune connections on the B-matrix only
[rPEB,BMR,BMA] = spm_dcm_bmr_all(PEB,{'B'});

% Check the only surviving B-matrix (group difference) parameter is R1->R2
b_reduced = BMA.Cp(:,2);
expected  = logical([0 0 0 0 1 0]);
testCase.assertTrue(all(b_reduced(~expected) < 1e-4));
testCase.assertTrue(all(b_reduced(expected) > 1e-4));

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');