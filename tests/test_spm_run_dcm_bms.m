function tests = test_spm_run_dcm_bms
% Unit Tests for config/spm_run_dcm_bms. Tests are provided with and
% without evidence for a particular model with artificially generated free
% energies. Additionally, tests are included using real DCM files for
% software testing.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_run_dcm_bms.m 7103 2017-06-08 14:39:32Z peter $

tests = functiontests(localfunctions);

function setup(testCase)
% Prepare output directory
out_dir = get_output_dir();

% Delete existing files if they exist
if exist(fullfile(out_dir,'BMS.mat'),'file')
    delete(fullfile(out_dir,'BMS.mat'));
end
if exist(fullfile(out_dir,'F.mat'),'file')
    delete(fullfile(out_dir,'F.mat'));
end

% -------------------------------------------------------------------------
function test_ffx_fmri_dcm(testCase)
% Tests that the BMS can run on an fMRI DCM (model 2 should win)

models_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region', 'models');

% Select DCMs (2 per subject)
P1 = spm_select('FPList',models_path,'DCM_.*_m1\.mat$');
P2 = spm_select('FPList',models_path,'DCM_.*_m2\.mat$');
P  = [cellstr(P1)  cellstr(P2)];

% Run FFX BMS
run_bms_dcm_files('FFX', P);

% Check
load(fullfile(get_output_dir(), 'BMS.mat'));
testCase.assertTrue(BMS.DCM.ffx.model.post(2) > 0.99);
% -------------------------------------------------------------------------
function test_ffx_csd_dcm(testCase)
% Tests that the BMS can run on a CSD DCM

models_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'MEEG');

% Select DCMs (2 per subject)
P1 = spm_select('FPList',models_path,'DCM_CSD.mat');
P2 = spm_select('FPList',models_path,'DCM_CSD.mat');
P  = [cellstr(P1)  cellstr(P2)];

% Run FFX BMS
run_bms_dcm_files('FFX', P);

% Check
load(fullfile(get_output_dir(), 'BMS.mat'));
testCase.assertTrue(BMS.DCM.ffx.model.post(2) == 0.5);
% -------------------------------------------------------------------------
function test_rfx_no_evidence(testCase)
% Tests RFX BMS in the context of equal model probabilities

% Setup
import matlab.unittest.constraints.*
rng('default');
rng(1);

n = 20;      % Subjects
models = 4;  % Models per subject

% Frequency of each model in the group
r = [0.25 0.25 0.25 0.25];

num_subjects_per_model = n .* r;

% Build subjects x models matrix of model assignments
m = [];
for model = 1:models  
    row = zeros(1,models); 
    row(model) = 1;
    
    m = [m; repmat(row, num_subjects_per_model(model), 1)];
end

% Free energy of the generative model in each subject
mean_F = -10; 
std_F  = 0.5;
F_gen  = mean_F + std_F .* randn(sum(m(:) > 0),1);

% Build free energy matrix (mean -13, STD 1)
F    = -13 + randn(n,4);
F(m > 0) = F_gen;
save('output/F.mat','F');

% Run
run_bms_Fmatrix('RFX');

load('output/BMS.mat');

% Check expected frequencies are within 10%
actual = BMS.DCM.rfx.model.exp_r;
testCase.verifyThat(actual, IsEqualTo(r, 'Within', AbsoluteTolerance(0.1) ) );

% Test BOR
actual = BMS.DCM.rfx.model.bor;
testCase.verifyThat(actual, IsGreaterThanOrEqualTo(0.9 ) );

% -------------------------------------------------------------------------
function test_rfx_strong_evidence(testCase)
% Tests RFX BMS in the context of an effect

% Setup
import matlab.unittest.constraints.*
rng('default');
rng(1);

n = 20;      % Subjects
models = 4;  % Models per subject

% Frequency of each model in the group
r = [0.60 0.20 0.15 0.05];

num_subjects_per_model = n .* r;

% Build subjects x models matrix of model assignments
m = [];
for model = 1:models  
    row = zeros(1,models); 
    row(model) = 1;
    
    m = [m; repmat(row, num_subjects_per_model(model), 1)];
end

% Free energy of the generative model in each subject
mean_F = -10; 
std_F  = 0.5;
F_gen  = mean_F + std_F .* randn(sum(m(:) > 0),1);

% Build free energy matrix (mean -13, STD 1)
F    = -13 + randn(n,4);
F(m > 0) = F_gen;
save('output/F.mat','F');

% Run
run_bms_Fmatrix('RFX');

load('output/BMS.mat');

% Check expected frequencies are within 10%
actual = BMS.DCM.rfx.model.exp_r;
testCase.verifyThat(actual, IsEqualTo(r, 'Within', AbsoluteTolerance(0.1) ) );

% Test BOR
actual = BMS.DCM.rfx.model.bor;
testCase.verifyThat(actual, IsLessThanOrEqualTo(0.1 ) );

% -------------------------------------------------------------------------
function test_ffx_strong_evidence(testCase)
% Tests FFX in the context of strong evidence

import matlab.unittest.constraints.*

% Free energies for model 1 and model 2
F = [-10.3 -10.3 -10.3 -10.3 -10.3;
     -9.7  -9.7  -9.7  -9.7  -9.7]';
save('output/F.mat','F');

% Run
run_bms_Fmatrix('FFX');

% Check    
load('output/BMS.mat');    
actual_group_log_bf   = BMS.DCM.ffx.SF - min(BMS.DCM.ffx.SF);
actual_PP             = BMS.DCM.ffx.model.post(2);

testCase.verifyThat(actual_group_log_bf(2), ...
    IsEqualTo(3, 'Within', AbsoluteTolerance(0.01) ) );

testCase.verifyThat(actual_PP, ...
    IsEqualTo(0.95, 'Within', AbsoluteTolerance(0.01)));

% -------------------------------------------------------------------------
function test_ffx_no_evidence(testCase)
% Tests FFX in the context of no evidence

import matlab.unittest.constraints.*

% Free energies for model 1 and model 2
F = [-10 -10 -10 -10 -10
     -10 -10 -10 -10 -10]';
save('output/F.mat','F');

% Run
run_bms_Fmatrix('FFX');

% Check    
load('output/BMS.mat');    
actual_group_log_bf   = BMS.DCM.ffx.SF - min(BMS.DCM.ffx.SF);
actual_PP             = BMS.DCM.ffx.model.post(2);

testCase.verifyThat(actual_group_log_bf(2), ...
    IsEqualTo(0, 'Within', AbsoluteTolerance(0.01) ) );

testCase.verifyThat(actual_PP, ...
    IsEqualTo(0.5, 'Within', AbsoluteTolerance(0.01)));

% -------------------------------------------------------------------------
function run_bms_Fmatrix(method)
% Run BMS using saved log evidence matrix
clear matlabbatch;
matlabbatch{1}.spm.dcm.bms.inference.dir = {'output'};
matlabbatch{1}.spm.dcm.bms.inference.sess_dcm = {};
matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
matlabbatch{1}.spm.dcm.bms.inference.load_f = {'output/F.mat'};
matlabbatch{1}.spm.dcm.bms.inference.method = method;
matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
spm_jobman('run',matlabbatch);
% -------------------------------------------------------------------------
function run_bms_dcm_files(method, P)
% Run BMS using DCM .mat files
%
% P - subjects x models cell array

clear matlabbatch;
matlabbatch{1}.spm.dcm.bms.inference.dir = cellstr(get_output_dir());
matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
matlabbatch{1}.spm.dcm.bms.inference.method = method;
matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
matlabbatch{1}.spm.dcm.bms.inference.verify_id = 0;

for s = 1:size(P,1)
    matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{s}.dcmmat = P(s,:)';
end
spm_jobman('run',matlabbatch);

% -------------------------------------------------------------------------
function out_dir = get_output_dir()
% Returns the directory for output files and creates it if needed
out_dir = fullfile( spm('Dir'), 'tests', 'output');
if ~exist(out_dir,'file')
    mkdir(out_dir);
end
