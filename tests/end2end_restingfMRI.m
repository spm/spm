function tests = end2end_restingfMRI
% End-to-end test for resting dataset
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% $Id: end2end_restingfMRI.m 7481 2018-11-09 15:36:57Z peter $

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function setupOnce(testCase)

% Prepare path
data_path = get_data_path();
rest_path  = fullfile(data_path,'spdcm');

if exist(rest_path,'file')
    % Clear out existing analysis files
    P = spm_select('FPList', fullfile(rest_path,'glm'),'.*');
    for i = 1:size(P,1)
        spm_unlink(P(i,:));
    end
    
    P = spm_select('FPList', fullfile(rest_path,'glm_corrected'),'.*');
    for i = 1:size(P,1)
        spm_unlink(P(i,:));
    end    
else
    mkdir(rest_path);
end

% Download dataset
spm_get_dataset('spm', 'spdcm', '', data_path);

% -------------------------------------------------------------------------
function test_CSD(testCase)

spm('defaults','fmri')

% Run GLM
run_glm();

% Specify DCMs
specify_dcms();

dcmmat = fullfile(get_data_path(),'spdcm','glm_corrected','DCM_full.mat');

% Estimate single DCM
clear matlabbatch;
matlabbatch{1}.spm.dcm.estimate.dcms.subj.dcmmat = cellstr(dcmmat);
matlabbatch{1}.spm.dcm.estimate.output.separate = struct([]);
matlabbatch{1}.spm.dcm.estimate.est_type = 3;
matlabbatch{1}.spm.dcm.estimate.fmri.analysis = 'csd';
spm_jobman('run',matlabbatch);

% Load created DCM
DCM = load(dcmmat);
DCM = DCM.DCM;

% Check model fit (full model)
DCM = spm_dcm_fmri_check(DCM,true);
exp_var = DCM.diagnostics(1);
max_A   = DCM.diagnostics(2);
testCase.assertTrue(exp_var > 95);
testCase.assertTrue(max_A > 1);

% -------------------------------------------------------------------------
function run_glm()
% Runs GLM and VOI extraction

data_path = fullfile(get_data_path(), 'spdcm');

% Initialise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline',1);

f = spm_select('FPList', fullfile(data_path,'func'), '^sw.*\.img$');


RT = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL GLM FOR EXTRACTING WM / CSF REGRESSORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

glmdir = fullfile(data_path,'glm');
if ~exist(glmdir,'file'), mkdir(glmdir); end

clear matlabbatch;

% SPM specification
matlabbatch{1}.spm.stats.fmri_spec.dir          = cellstr(glmdir);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT    = RT;
matlabbatch{1}.spm.stats.fmri_spec.sess.scans   = cellstr(f);

% SPM estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(glmdir,'SPM.mat'));

% ROI extraction
matlabbatch{3}.spm.util.voi.spmmat  = cellstr(fullfile(glmdir,'SPM.mat'));
matlabbatch{3}.spm.util.voi.adjust  = NaN;
matlabbatch{3}.spm.util.voi.session = 1;
matlabbatch{3}.spm.util.voi.name    = 'CSF';
matlabbatch{3}.spm.util.voi.roi{1}.sphere.centre     = [ 0 -40 -5];
matlabbatch{3}.spm.util.voi.roi{1}.sphere.radius     = 6;
matlabbatch{3}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
matlabbatch{3}.spm.util.voi.roi{2}.mask.image        = cellstr(fullfile(glmdir,'mask.nii'));
matlabbatch{3}.spm.util.voi.expression = 'i1 & i2';

matlabbatch{4} = matlabbatch{3};
matlabbatch{4}.spm.util.voi.name = 'WM';
matlabbatch{4}.spm.util.voi.roi{1}.sphere.centre = [0 -24 -33]; 

spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECOND GLM INCLUDING WM / CSF REGRESSORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

glmdir = fullfile(data_path,'glm_corrected');
if ~exist(glmdir,'file'), mkdir(glmdir); end

clear matlabbatch;

% SPM specification
matlabbatch{1}.spm.stats.fmri_spec.dir          = cellstr(glmdir);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT    = RT;
matlabbatch{1}.spm.stats.fmri_spec.sess.scans   = cellstr(f);
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {
    fullfile(data_path,'glm','rp_rest0000.txt'),...
    fullfile(data_path,'glm','VOI_CSF_1.mat'),...
    fullfile(data_path,'glm','VOI_WM_1.mat'),...
    }';

% SPM estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(glmdir,'SPM.mat'));

% ROI extraction
matlabbatch{3}.spm.util.voi.spmmat  = cellstr(fullfile(glmdir,'SPM.mat'));
matlabbatch{3}.spm.util.voi.adjust  = NaN;
matlabbatch{3}.spm.util.voi.session = 1;
matlabbatch{3}.spm.util.voi.name    = 'PCC';
matlabbatch{3}.spm.util.voi.roi{1}.sphere.centre = [0 -52 26];
matlabbatch{3}.spm.util.voi.roi{1}.sphere.radius = 8;
matlabbatch{3}.spm.util.voi.roi{2}.mask.image    = cellstr(fullfile(glmdir,'mask.nii'));
matlabbatch{3}.spm.util.voi.expression = 'i1 & i2';

matlabbatch{4} = matlabbatch{3};
matlabbatch{4}.spm.util.voi.name = 'mPFC';
matlabbatch{4}.spm.util.voi.roi{1}.sphere.centre = [3 54 -2];

matlabbatch{5} = matlabbatch{3};
matlabbatch{5}.spm.util.voi.name = 'LIPC';
matlabbatch{5}.spm.util.voi.roi{1}.sphere.centre = [-50 -63 32];

matlabbatch{6} = matlabbatch{3};
matlabbatch{6}.spm.util.voi.name = 'RIPC';
matlabbatch{6}.spm.util.voi.roi{1}.sphere.centre = [48 -69 35];

spm_jobman('run',matlabbatch);

% -------------------------------------------------------------------------
function specify_dcms()

% SPM
glm_dir = fullfile(get_data_path(),'spdcm','glm_corrected');
SPM = fullfile(glm_dir,'SPM.mat');

% VOIs
xY  = {fullfile(glm_dir,'VOI_PCC_1.mat');
       fullfile(glm_dir,'VOI_mPFC_1.mat');
       fullfile(glm_dir,'VOI_LIPC_1.mat');
       fullfile(glm_dir,'VOI_RIPC_1.mat')};

n   = 4; % Num regions
nu  = 1; % Num inputs

% Connectivity matrices
a  = ones(n,n);
b  = zeros(n,n,nu);
c  = zeros(n,nu);
d  = zeros(n,n,0);

% Specify model 1 (full)
s = struct();
s.name       = 'full';
s.u          = [1 1 1]';
s.delays     = repmat(2,1,n);
s.TE         = 0.04;
s.nonlinear  = false;
s.two_state  = false;
s.stochastic = false;
s.centre     = false;
s.induced    = 1;
s.a          = a;
s.b          = b;
s.c          = c;
s.d          = d;
spm_dcm_specify(SPM,xY,s);

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'output');