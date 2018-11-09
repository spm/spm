function tests = test_spm_cfg_dcm_est
% Unit Tests for test_spm_cfg_dcm_est (DCM model estimation batch)
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_cfg_dcm_est.m 7479 2018-11-09 14:17:33Z peter $

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function setup(testCase)
% Delete output directory before each test
[input_path, path] = prepare_paths(0);

% Paths to delete
paths = {fullfile(path,'tmp'); path};

for p = 1:length(paths)
    if exist(paths{p},'file')
        delete(fullfile(paths{p},'*.mat'));
        rmdir(paths{p});
    end
end

% -------------------------------------------------------------------------
function test_select_by_gcm(testCase)
% Test selecting a group DCM file (GCM)

% Prepare paths
[input_path, path] = prepare_paths();

% Expected output
expected = fullfile(path,'GCM_test.mat');

% Input
gcm_file = fullfile(input_path,'GCM_simulated.mat'); 

% Complete & run batch
matlabbatch{1}.spm.dcm.estimate.dcms.gcmmat = cellstr(gcm_file);
matlabbatch{1}.spm.dcm.estimate.output.single.dir = cellstr(path);
matlabbatch{1}.spm.dcm.estimate.output.single.name = 'test';
matlabbatch{1}.spm.dcm.estimate.est_type = 4;
matlabbatch{1}.spm.dcm.estimate.fmri.analysis = 'time';
out = spm_jobman('run',matlabbatch);

% Check output created
actual = out{1}.gcmmat{1};
testCase.assertEqual(actual,expected);
testCase.assertTrue(exist(expected,'file') > 0);

% Check output contents is correct
assert_gcms_match(actual,gcm_file,testCase);

% -------------------------------------------------------------------------
function test_select_by_model(testCase)
% Test selecting DCMs model-by-model

% Prepare paths
[input_path, path] = prepare_paths();

% Expected output
expected = fullfile(path,'GCM_test.mat');

nm = 4;
ns = 30;

% Gather models
for m = 1:nm    
    P = {};
    for s = 1:ns    
        str = sprintf('DCM_s%d_m%d.mat',s,m);
        P   = [P; cellstr(fullfile(input_path,str))];            
    end
    matlabbatch{1}.spm.dcm.estimate.dcms.model(m).dcmmat = P;
end

% Complete & run batch
matlabbatch{1}.spm.dcm.estimate.output.single.dir = cellstr(path);
matlabbatch{1}.spm.dcm.estimate.output.single.name = 'test';
matlabbatch{1}.spm.dcm.estimate.est_type = 4;
matlabbatch{1}.spm.dcm.estimate.fmri.analysis = 'time';
out = spm_jobman('run',matlabbatch);

% Check output created
actual = out{1}.gcmmat{1};
testCase.assertEqual(actual,expected);
testCase.assertTrue(exist(expected,'file') > 0);

% Check output contents is correct
gcm_file = fullfile(input_path,'GCM_simulated.mat'); 
assert_gcms_match(actual,gcm_file,testCase);

% -------------------------------------------------------------------------
function test_select_by_subject(testCase)
% Test selecting DCMs subject-by-subject

% Prepare paths
[input_path, path] = prepare_paths();

% Expected output
expected = fullfile(path,'GCM_test.mat');

ns = 30;

% Gather models
for s = 1:ns    
    str = sprintf('^DCM_s%d_.*mat$',s);
    P   = spm_select('FPList',input_path,str);
    
    matlabbatch{1}.spm.dcm.estimate.dcms.subj(s).dcmmat = cellstr(P);
end

% Complete & run batch
matlabbatch{1}.spm.dcm.estimate.output.single.dir = cellstr(path);
matlabbatch{1}.spm.dcm.estimate.output.single.name = 'test';
matlabbatch{1}.spm.dcm.estimate.est_type = 4;
matlabbatch{1}.spm.dcm.estimate.fmri.analysis = 'time';
out = spm_jobman('run',matlabbatch);

% Check output created
actual = out{1}.gcmmat{1};
testCase.assertEqual(actual,expected);
testCase.assertTrue(exist(expected,'file') > 0);

% Check output contents is correct
gcm_file = fullfile(input_path,'GCM_simulated.mat'); 
assert_gcms_match(actual,gcm_file,testCase);

% -------------------------------------------------------------------------
function test_separate_dcm_output(testCase)
% Test that outputing separate DCMs gives the correct results

% Prepare paths
[input_path, path] = prepare_paths();

% Copy all DCMs into a temporary folder
tmp_path = fullfile(path, 'tmp');
mkdir(tmp_path);
copyfile( fullfile(input_path,'DCM_*.mat'), tmp_path); 

ns = 30;

% Gather models
P = {};
for s = 1:ns    
    str    = sprintf('^DCM_s%d_.*mat$',s);
    P(s,:) = cellstr(spm_select('FPList',tmp_path,str));
    
    matlabbatch{1}.spm.dcm.estimate.dcms.subj(s).dcmmat = P(s,:)';
end

% Complete & run batch
matlabbatch{1}.spm.dcm.estimate.output.separate = struct([]);
matlabbatch{1}.spm.dcm.estimate.est_type = 4;
matlabbatch{1}.spm.dcm.estimate.fmri.analysis = 'time';
out = spm_jobman('run',matlabbatch);

% Check output created
actual   = out{1}.dcmmat;
expected = P(:,1);
testCase.assertEqual(actual,expected);

% Check output contents is correct
gcm_file = fullfile(input_path,'GCM_simulated.mat'); 
GCM = spm_dcm_load(gcm_file);
assert_gcms_match(actual,GCM(:,1),testCase);

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');

% -------------------------------------------------------------------------
function [input_path, path] = prepare_paths(do_create)
% Prepare input and output paths. If do_create, create out directory

if nargin < 1
    do_create = 1;
end

data_path   = get_data_path();
input_path  = fullfile(data_path,'models');
path = fullfile(data_path,'out');

if do_create && ~exist(path,'file')
    mkdir(path);
end
% -------------------------------------------------------------------------
function assert_gcms_match(actual,expected,testCase)
% Check that two GCM cell arrays match in terms of contents

% Load
if ischar(actual)
    actual = load(actual);
    actual = actual.GCM;
end

if ischar(expected)
    expected = load(expected);
    expected = expected.GCM;
end

% Convert GCM of filenames -> GCM of DCMs if needed
if ischar(actual{1}),   actual   = spm_dcm_load(actual); end
if ischar(expected{1}), expected = spm_dcm_load(expected); end

% Check sizes match
testCase.assertEqual(size(actual),size(expected));

% Check each DCM matches
nm = numel(actual);

for m = 1:nm
    actual_id   = spm_data_id(actual{m});
    expected_id = spm_data_id(expected{m});
    testCase.assertEqual(actual_id,expected_id,...
        'AbsTol',0.001,sprintf('GCM index %d doesn''t match',m));
end