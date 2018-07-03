function tests = test_spm_dcm_fit
% Unit Tests for spm_dcm_fit
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dcm_fit.m 7364 2018-07-03 14:02:46Z peter $

tests = functiontests(localfunctions);

function test_fmri_noparfor(testCase)

data_path = get_data_path();

% Load DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;
GCM = GCM(1:2,1);

% Remove old posteriors
for i = 1:length(GCM)
    GCM{i} = rmfield(GCM{i},'Ep');
    GCM{i} = rmfield(GCM{i},'F');
end

% Estimate
GCM = spm_dcm_fit(GCM);

% Check posteriors
for i = 1:length(GCM)
    testCase.assertTrue(isfield(GCM{i},'Ep'));
    testCase.assertTrue(isfield(GCM{i},'F'));
end

% -------------------------------------------------------------------------
function test_fmri_parfor(testCase)

data_path = get_data_path();

% Load DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;
GCM = GCM(1:2,1);

% Remove old posteriors
for i = 1:length(GCM)
    GCM{i} = rmfield(GCM{i},'Ep');
    GCM{i} = rmfield(GCM{i},'F');
end

% Estimate
GCM = spm_dcm_fit(GCM, true);

% Check posteriors
for i = 1:length(GCM)
    testCase.assertTrue(isfield(GCM{i},'Ep'));
    testCase.assertTrue(isfield(GCM{i},'F'));
end

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');