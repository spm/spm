function tests = test_spm_dcm_fmri_check
% Unit Tests for spm_dcm_fmri_check
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dcm_fmri_check.m 6788 2016-04-28 11:45:23Z peter $

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_DCM(testCase)
% Simply tests the function doesn't crash with one DCM

data_path = get_data_path();
dcm_file  = fullfile(data_path,'DCM_s1_m1.mat');

DCM = spm_dcm_fmri_check(dcm_file, true);

assert(isfield(DCM,'diagnostics'));
assert(numel(DCM.diagnostics) == 3);

% -------------------------------------------------------------------------
function test_GCM(testCase)
% Simply tests the function doesn't crash with a GCM array

data_path = get_data_path();
gcm_file  = fullfile(data_path,'GCM_simulated.mat');

GCM = spm_dcm_fmri_check(gcm_file, true);

assert(iscell(GCM));
assert(isfield(GCM{1},'diagnostics'));
assert(numel(GCM{1}.diagnostics) == 3);

% -------------------------------------------------------------------------
function data_path = get_data_path()
data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region', 'models');

