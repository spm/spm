function tests = test_spm_dcm_peb_review
% Unit Tests for test_spm_dcm_peb_review. Simply ensures that the GUI
% doesn't crash with different inputs.
%__________________________________________________________________________

% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_with_peb(testCase)

data_path = get_data_path();

PEB = fullfile(data_path,'PEB_test.mat');
DCM = fullfile(data_path,'models','DCM_s1_m1.mat');

spm_dcm_peb_review(PEB,DCM);
close all;
% -------------------------------------------------------------------------
function test_with_unestimated_dcm(testCase)

data_path = get_data_path();

PEB = fullfile(data_path,'PEB_test.mat');
DCM = fullfile(data_path,'models','DCM_s1_m1.mat');

DCM = load(DCM);
DCM = DCM.DCM;

% Create simpler fMRI DCM with minimal structures
DCM2 = struct();
DCM2.a = DCM.a;
DCM2.b = DCM.b;
DCM2.c = DCM.c;
DCM2.d = DCM.d;
DCM2.xY = DCM.xY;
DCM2.Y  = DCM.Y;
DCM2.U  = DCM.U;
DCM2.options=DCM.options;

spm_dcm_peb_review(PEB,DCM2);
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
