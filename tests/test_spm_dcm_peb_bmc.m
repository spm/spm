function tests = test_spm_dcm_peb_bmc
% Unit Tests for test_spm_dcm_peb_bmc
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dcm_peb_bmc.m 7327 2018-06-07 13:30:01Z peter $

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_model_search(testCase)

data_path = get_data_path();

PEB = load(fullfile(data_path,'PEB_test.mat'));
PEB = PEB.PEB;

[BMA,BMR] = spm_dcm_peb_bmc(PEB(1));

% Get BMA parameters
[n,m] = size(PEB(1).Ep);
Ep = reshape(full(BMA.Ep),[n,m]);
Pp = reshape(full(BMA.Pp),[n,m]);

effect     = 2; % group difference
connection = 5; % r1->r2

% There should be an effect of group on the connection from R1->R2 of 0.2Hz 
testCase.assertEqual(Ep(connection,effect),0.2, 'AbsTol', 0.05);
testCase.assertTrue(Pp(connection,effect) > 0.95);

% There should be no effect of group elsewhere
Pp_others = Pp; 
Pp_others(connection,:) = [];
testCase.assertTrue(all(Pp_others(connection,effect) < 0.95));

close all;

% -------------------------------------------------------------------------
function test_specific_model_comparison(testCase)

data_path = get_data_path();

% Get PEB full model
PEB = load(fullfile(data_path,'PEB_test.mat'));
PEB = PEB.PEB(1);

% Get template models
GCM_templates = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM_templates = GCM_templates.GCM;

% Run model comparison
[BMA,BMR] = spm_dcm_peb_bmc(PEB(1),GCM_templates(1,:));

% Get BMA parameters
[n,m] = size(PEB(1).Ep);
Ep = reshape(full(BMA.Ep),[n,m]);

effect     = 2; % group difference
connection = 5; % r1->r2

% There should be common effects on the forward connection only
testCase.assertTrue(BMA.Pw(1) > 0.95);
testCase.assertTrue(BMA.Pw(2) < 0.95);

% There should be an effect of group on the connection from R1->R2 of 0.2Hz 
testCase.assertEqual(Ep(connection,effect),0.2, 'AbsTol', 0.05);
testCase.assertTrue(BMA.Px(1) > 0.95);
testCase.assertTrue(BMA.Px(2) < 0.95);

close all;

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');