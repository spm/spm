function tests = test_regress_spm_dcm_fmri
% Regression test for DCM for fMRI including timseries extraction
% % This script analyses the Attention to Visual Motion fMRI dataset
% available from the SPM website using DCM:
%   http://www.fil.ion.ucl.ac.uk/spm/data/attention/
% as described in the SPM docs website:
%   https://www.fil.ion.ucl.ac.uk/spm/docs/tutorials/dcm/dcm_fmri_first_level_gui/
%__________________________________________________________________________

% Copyright (C) 2024 Wellcome Centre for Human Neuroimaging

tests = functiontests(localfunctions);


function test_regress_attention_dcm(testCase)

data_path = fullfile(spm('Dir'),'tests','data','attention');

% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM SPECIFICATION, ESTIMATION & INFERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factors = load(fullfile(data_path,'factors.mat'));

f = spm_select('FPList', fullfile(data_path,'functional'), '^snf.*\.img$');

clear matlabbatch

% OUTPUT DIRECTORY
%--------------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'GLM';

% MODEL SPECIFICATION
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(fullfile(data_path,'GLM'));
matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{2}.spm.stats.fmri_spec.timing.RT    = 3.22;
matlabbatch{2}.spm.stats.fmri_spec.sess.scans            = cellstr(f);
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).name     = 'Photic';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).onset    = [factors.phot];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).duration = 10;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).name     = 'Motion';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).onset    = [factors.mot];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).duration = 10;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).name     = 'Attention';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).onset    = [factors.att];
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).duration = 10;

% MODEL ESTIMATION
%--------------------------------------------------------------------------
matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));

% INFERENCE
%--------------------------------------------------------------------------
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VOLUMES OF INTEREST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

% EXTRACTING TIME SERIES: V5
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DYNAMIC CAUSAL MODELLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select SPM.mat
spm_file = fullfile(data_path,'GLM','SPM.mat');

% Select VOIs
xY = cell(3,1);
xY{1} = fullfile(data_path,'GLM','VOI_V1_1.mat');
xY{2} = fullfile(data_path,'GLM','VOI_V5_1.mat');
xY{3} = fullfile(data_path,'GLM','VOI_SPC_1.mat');

% Indices of the brain regions in the model
V1  = 1; 
V5  = 2;
SPC = 3;

% DCM settings
n   = 3;    % number of regions
nu  = 3;    % number of inputs (experimental conditions)
TR  = 3.22; % volume repetition time (seconds)
TE  = 0.04; % echo time (seconds)

% Experimental conditions to include from the SPM.
cond = struct();
cond(1).name = 'Photic';
cond(2).name = 'Motion';
cond(3).name = 'Attention';

% Indices of the conditions in the model
PHOTIC = 1;
MOTION = 2;
ATTENTION = 3;

% Connectivity matrices
a  = ones(n,n);
b  = zeros(n,n,nu);
c  = zeros(n,nu);
d  = zeros(n,n,0);

% A-matrix connectivity
a(V1,SPC) = 0; % SPC->V1
a(SPC,V1) = 0; % V1->SPC

% Modulatory input: motion on V1->V5
b(V5,V1,MOTION) = 1;

% Driving input: photic
c(V1,PHOTIC) = 1;

% DCM options
s = struct();
s.cond       = cond;
s.delays     = repmat(TR/2, 1, n);
s.TE         = TE;
s.nonlinear  = false;
s.two_state  = false;
s.stochastic = false;
s.centre     = true;
s.induced    = 0;
s.a          = a;
s.b          = b;
s.c          = c;
s.d          = d;

% Cell array to store DCMs [1 subject x 2 models]
GCM = cell(1,2);

% Specify model 1: forward model, with attention on V1->V5
s.name = 'mod_fwd';
s.b(V5,V1,ATTENTION) = 1;
GCM{1,1} = spm_dcm_specify(spm_file,xY,s);

% Specify model 2: backward model, with attention on SPC->V5
s.name = 'mod_bwd';
s.b = b;
s.b(V5,SPC,ATTENTION);
GCM{1,2} = spm_dcm_specify(spm_file,xY,s);

% Suppress output
spm_get_defaults('dcm.verbose',false);

% Estimate the models
GCM = spm_dcm_fit(GCM);

% Bayesian model comparison
post = spm_dcm_bmc(GCM);

% Get the explained variance
DCM = spm_dcm_fmri_check(GCM{1},true);
exp_var = DCM.diagnostics(1);
max_connection = DCM.diagnostics(2);

% Check the forward model was the winner (Pp > 95%)
testCase.verifyTrue(post(1) > 0.95); 

% Check the explained variance was high (> 75%)
testCase.verifyTrue(exp_var > 75);

% Check the largest neural connection was non-trivial (> 0.5Hz)
testCase.verifyTrue(max_connection > 0.5);

