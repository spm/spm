function tests = end2end_attention
% End-to-end test for attention dataset
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% $Id: end2end_attention.m 7264 2018-02-22 14:43:47Z peter $

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function setupOnce(testCase)

% Prepare path
data_path = get_data_path();
att_path  = fullfile(data_path,'attention');

if exist(att_path,'file')
    % Clear out existing analysis files
    P = spm_select('FPList', fullfile(att_path,'GLM'),'.*');
    for i = 1:size(P,1)
        spm_unlink(P(i,:));
    end
else
    mkdir(att_path);
end

% Download dataset
spm_get_dataset('spm', 'attention', '', data_path);

% -------------------------------------------------------------------------
function test_attention(testCase)

% Run GLM
run_glm();

% Specify DCMs
specify_dcms();

% Select DCMs
glm_dir = fullfile(get_data_path(), 'attention', 'GLM');
P = {fullfile(glm_dir,'DCM_full.mat');
     fullfile(glm_dir,'DCM_fwd.mat');
     fullfile(glm_dir,'DCM_bwd.mat')};

% Estimate DCMs (VL + BMR)
clear matlabbatch;
matlabbatch{1}.spm.dcm.estimate.dcms.subj.dcmmat = P;
matlabbatch{1}.spm.dcm.estimate.output.single.dir = cellstr(glm_dir);
matlabbatch{1}.spm.dcm.estimate.output.single.name = 'full_fwd_bwd';
matlabbatch{1}.spm.dcm.estimate.est_type = 1;
matlabbatch{1}.spm.dcm.estimate.fmri.analysis = 'time';
spm_jobman('run',matlabbatch);

% Load GCM
gcm_file = fullfile(glm_dir,'GCM_full_fwd_bwd.mat');
GCM = spm_dcm_load(gcm_file);

% Compare models
post = spm_dcm_bmc(GCM);

% Check model 2 was the best
[maxP,idx] = max(post);
testCase.assertTrue(idx == 2);
testCase.assertTrue(maxP > 0.9);

% Check model fit (full model)
GCM = spm_dcm_fmri_check(GCM,true);
DCM = GCM{1};
exp_var = DCM.diagnostics(1);
max_A   = DCM.diagnostics(2);
testCase.assertTrue(exp_var > 85);
testCase.assertTrue(max_A > 0.4);

% Check the BOLD impulse response kernel looks sensible
x = (1:DCM.M.N)*DCM.M.dt;
y = DCM.H1(:,1,1);
[max_y, idx] = max(y);
max_x = x(idx);
testCase.assertTrue(max_x > 5 & max_x < 7);
testCase.assertTrue(max_y > 0);

% -------------------------------------------------------------------------
function run_glm()
% Runs GLM and VOI extraction

data_path = fullfile(get_data_path(), 'attention');

% Initialise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline',1);

factors = load(fullfile(data_path,'factors.mat'));
f = spm_select('FPList', fullfile(data_path,'functional'), '^snf.*\.img$');

clear matlabbatch

% OUTPUT DIRECTORY
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'GLM';

% MODEL SPECIFICATION
matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(fullfile(data_path,'GLM'));
matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{2}.spm.stats.fmri_spec.timing.RT    = 3.22;
matlabbatch{2}.spm.stats.fmri_spec.sess.scans            = cellstr(f);
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).name     = 'Photic';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).onset    = [factors.att factors.natt factors.stat];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).duration = 10;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).name     = 'Motion';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).onset    = [factors.att factors.natt];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).duration = 10;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).name     = 'Attention';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).onset    = [factors.att];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).duration = 10;

% MODEL ESTIMATION
matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));

% INFERENCE
matlabbatch{4}.spm.stats.con.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{4}.spm.stats.con.consess{1}.fcon.name = 'Effects of Interest';
matlabbatch{4}.spm.stats.con.consess{1}.fcon.weights = eye(3);
matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'Photic';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [1 0 0];
matlabbatch{4}.spm.stats.con.consess{3}.tcon.name = 'Motion';
matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [0 1 0];
matlabbatch{4}.spm.stats.con.consess{4}.tcon.name = 'Attention';
matlabbatch{4}.spm.stats.con.consess{4}.tcon.weights = [0 0 1];

spm_jobman('run',matlabbatch);

% VOLUMES OF INTEREST
clear matlabbatch

% EXTRACTING TIME SERIES: V5
matlabbatch{1}.spm.util.voi.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{1}.spm.util.voi.adjust = 1;  % "effects of interest" F-contrast
matlabbatch{1}.spm.util.voi.session = 1; % session 1
matlabbatch{1}.spm.util.voi.name = 'V5';
matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {''}; % using SPM.mat above
matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 3;  % "Motion" T-contrast
matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'FWE';
matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.05;
matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask.contrast = 4; % "Attention" T-contrast
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask.thresh = 0.05;
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask.mtype = 0; % inclusive
matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = [-36 -87 -3];
matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 8;
matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';

% EXTRACTING TIME SERIES: V1
matlabbatch{2}.spm.util.voi.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{2}.spm.util.voi.adjust = 1;  % "effects of interest" F-contrast
matlabbatch{2}.spm.util.voi.session = 1; % session 1
matlabbatch{2}.spm.util.voi.name = 'V1';
matlabbatch{2}.spm.util.voi.roi{1}.spm.spmmat = {''}; % using SPM.mat above
matlabbatch{2}.spm.util.voi.roi{1}.spm.contrast = 2;  % "Photic" T-contrast
matlabbatch{2}.spm.util.voi.roi{1}.spm.threshdesc = 'FWE';
matlabbatch{2}.spm.util.voi.roi{1}.spm.thresh = 0.05;
matlabbatch{2}.spm.util.voi.roi{1}.spm.extent = 0;
matlabbatch{2}.spm.util.voi.roi{2}.sphere.centre = [0 -93 18];
matlabbatch{2}.spm.util.voi.roi{2}.sphere.radius = 8;
matlabbatch{2}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{2}.spm.util.voi.expression = 'i1 & i2';

% EXTRACTING TIME SERIES: SPC
matlabbatch{3}.spm.util.voi.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{3}.spm.util.voi.adjust = 1;  % "effects of interest" F-contrast
matlabbatch{3}.spm.util.voi.session = 1; % session 1
matlabbatch{3}.spm.util.voi.name = 'SPC';
matlabbatch{3}.spm.util.voi.roi{1}.spm.spmmat = {''}; % using SPM.mat above
matlabbatch{3}.spm.util.voi.roi{1}.spm.contrast = 4;  % "Attention" T-contrast
matlabbatch{3}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
matlabbatch{3}.spm.util.voi.roi{1}.spm.thresh = 0.001;
matlabbatch{3}.spm.util.voi.roi{1}.spm.extent = 0;
matlabbatch{3}.spm.util.voi.roi{2}.sphere.centre = [-27 -84 36];
matlabbatch{3}.spm.util.voi.roi{2}.sphere.radius = 8;
matlabbatch{3}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{3}.spm.util.voi.expression = 'i1 & i2';

spm_jobman('run',matlabbatch);

% -------------------------------------------------------------------------
function specify_dcms()

% SPM
glm_dir = fullfile(get_data_path(),'attention','GLM');
SPM = fullfile(glm_dir,'SPM.mat');

% VOIs
xY  = {fullfile(glm_dir,'VOI_V1_1.mat');
       fullfile(glm_dir,'VOI_V5_1.mat');
       fullfile(glm_dir,'VOI_SPC_1.mat')};

n   = 3; % Num regions
nu  = 3; % Num inputs

% Indices of regions & conditions
V1        = 1; 
V5        = 2; 
SPC       = 3;
PHOTIC    = 1; 
MOTION    = 2; 
ATTENTION = 3;

% Connectivity matrices
a = ones(n,n);
b = zeros(n,n,nu);
c = zeros(n,nu);
d = zeros(n,n,0);
a(V1,SPC)          = 0;
a(SPC,V1)          = 0;
b(V5,V1,MOTION)    = 1;
b(V5,V1,ATTENTION) = 1;
b(V1,V5,ATTENTION) = 1;
c(V1,PHOTIC)       = 1;

% Specify model 1 (full)
s = struct();
s.name       = 'full';
s.u          = [1 1 1]';
s.delays     = [3.22 3.22 3.22];
s.TE         = 0.04;
s.nonlinear  = false;
s.two_state  = false;
s.stochastic = false;
s.centre     = false;
s.induced    = 0;
s.a          = a;
s.b          = b;
s.c          = c;
s.d          = d;
spm_dcm_specify(SPM,xY,s);

% Specify model 2 (fwd)
b(V5,V1,ATTENTION) = 1;
b(V1,V5,ATTENTION) = 0;
s.name = 'fwd';
s.b    = b;
spm_dcm_specify(SPM,xY,s);

% Specify model 2 (bwd)
b(V5,V1,ATTENTION) = 0;
b(V1,V5,ATTENTION) = 1;
s.name = 'bwd';
s.b    = b;
spm_dcm_specify(SPM,xY,s);

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'output');