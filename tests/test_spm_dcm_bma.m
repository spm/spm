function tests = test_spm_dcm_bma
% Unit Tests for spm_dcm_bma
%__________________________________________________________________________

% Copyright (C) 2016-2024 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_bma(testCase)

data_path = get_data_path();

% Load DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;

% Just keep the first group of subjects
% (see /data/fMRI/simulated_2region/generate_2regions.m)
GCM = GCM(1:15,:);

% Run BMA
BMA = spm_dcm_bma(GCM);

% Check model probabilities: model 2 should win
testCase.assertTrue(BMA.P(2) > 0.95);

% Indices in the DCM
INPUT_MODULATOR = 2;
R1 = 1; % Region 1
R2 = 2; % Region 2

% Get posterior expected value and variance for modulatory connection
ep_actual = BMA.Ep.B(R2,R1,INPUT_MODULATOR);
vp_actual = BMA.Cp.B(R2,R1,INPUT_MODULATOR);

% Get arithmetic mean of the modulatory parameter in model 2
Ep = cellfun(@(x)x.Ep.B(R2,R1,INPUT_MODULATOR),GCM(:,2));
expected = mean(Ep);

% Compute probability of difference from arithmetic mean
Pp = 1 - spm_Ncdf(expected,abs(ep_actual),vp_actual);
testCase.assertTrue(Pp < 0.7);
% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');
