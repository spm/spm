function tests = test_spm_dcm_post_hoc
% Unit Tests for spm_dcm_post_hoc
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dcm_post_hoc.m 6492 2015-06-26 14:27:40Z guillaume $

tests = functiontests(localfunctions);


%--------------------------------------------------------------------------
function data_path = get_data_path()
data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'test_spm_dcm_post_hoc');


%--------------------------------------------------------------------------
function setup(testCase)
data_path = get_data_path();

% Expected outputs
artefacts = {fullfile(data_path,'DCM_BPA.mat');
             fullfile(data_path,'DCM_opt_fwd_bwd_simulated.mat')};
    
% Delete if exist
for i = 1:length(artefacts)
    spm_unlink(artefacts{i});   
end

% Initialize SPM
spm('defaults','fmri');
spm_get_defaults('cmdline',true);


%--------------------------------------------------------------------------
function test_on_simulated_attention_data(testCase)
import matlab.unittest.constraints.*
 
data_path = get_data_path();

% Load model matching the generative model
DCM_fwd = load(fullfile(data_path, 'DCM_fwd_simulated.mat'));
DCM_fwd = DCM_fwd.DCM;

% Load model with extra connection
DCM_fwd_bwd = load(fullfile(data_path, 'DCM_fwd_bwd_simulated.mat'));
DCM_fwd_bwd = DCM_fwd_bwd.DCM;

% Run post-hoc
DCM_BPA = spm_dcm_post_hoc(fullfile(data_path, 'DCM_fwd_bwd_simulated.mat'));
DCM_opt = load(fullfile(data_path,'DCM_opt_fwd_bwd_simulated.mat'));
DCM_opt = DCM_opt.DCM;

% Check artefacts were created
testCase.assertEqual( exist(fullfile(data_path,'DCM_BPA.mat'),'file'), 2 );
testCase.assertEqual( exist(fullfile(data_path,'DCM_opt_fwd_bwd_simulated.mat'),'file'), 2 );

% Check reduced priors match the generative model
pC_fwd = full(diag(DCM_fwd.M.pC));
pC_opt = full(diag(DCM_opt.M.pC));
testCase.assertTrue(all(pC_fwd == pC_opt));

% % Plot for visual inspection
% spm_figure('Getwin','BMC'); clf
% 
% subplot(2,2,1); spm_plot_ci(DCM_fwd.Ep,DCM_fwd.Cp);
% title('Fwd (generative) model'); xlabel('Parameter')
% 
% subplot(2,2,2); spm_plot_ci(DCM_fwd_bwd.Ep,DCM_fwd_bwd.Cp);
% title('Bwd (over-complex) model'); xlabel('Parameter')
% 
% subplot(2,2,3); spm_plot_ci(DCM_BPA.Ep,DCM_BPA.Cp);
% title('BPA'); xlabel('Parameter')
% 
% subplot(2,2,4); spm_plot_ci(DCM_opt.Ep,DCM_opt.Cp);
% title('Optimal model for subject'); xlabel('Parameter');