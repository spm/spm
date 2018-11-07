function tests = test_spm_dcm_peb_review
% Unit Tests for test_spm_dcm_peb_review. Simply ensures that the GUI
% doesn't crash with different inputs.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dcm_peb_review.m 7474 2018-11-07 12:45:06Z peter $

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_with_peb(testCase)

data_path = get_data_path();

PEB = fullfile(data_path,'PEB_test.mat');
DCM = fullfile(data_path,'models','DCM_s1_m1.mat');

spm_dcm_peb_review(PEB,DCM);
close all;
% -------------------------------------------------------------------------
function test_with_bma_search(testCase)

data_path = get_data_path();

BMA = fullfile(data_path,'BMA_search.mat');
DCM = fullfile(data_path,'models','DCM_s1_m1.mat');

spm_dcm_peb_review(BMA,DCM);
close all;
% -------------------------------------------------------------------------
function test_with_bma_specific_models(testCase)

data_path = get_data_path();

BMA = fullfile(data_path,'BMA_specific_models.mat');
DCM = fullfile(data_path,'models','DCM_s1_m1.mat');

spm_dcm_peb_review(BMA,DCM);
close all;
% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');