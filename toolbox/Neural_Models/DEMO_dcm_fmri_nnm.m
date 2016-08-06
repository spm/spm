function DCM = DEMO_dcm_fmri_nnm
% test routine to check current implementations of DCM for fMRI
%==========================================================================
%
% Options
%--------------------------------------------------------------------------
% DCM.options.two_state              % two regional populations (E and I)
% DCM.options.stochastic             % fluctuations on hidden states
% DCM.options.nonlinear              % interactions among hidden states
% DCM.options.nograph                % graphical display
% DCM.options.centre                 % mean-centre inputs
% DCM.options.P                      % starting estimates for parameters
% DCM.options.hidden                 % indices of hidden regions

% $Id: DEMO_dcm_fmri_nnm.m 6855 2016-08-06 10:06:35Z karl $

% tests of spatial models: 'ECD', 'LFP' or 'IMG'
%==========================================================================
try
    cd(fullfile(spm('Dir'),'tests','data','fMRI'))
catch
    cd('C:\home\spm\DCM\DCM fMRI')
end
load DCM_attention

DCM = rmfield(DCM,{'a','b','c','d','options'});

% spatial models
%==========================================================================
DCM.options.nmm        = 'CMC';    % two regional populations (E and I)
DCM.options.stochastic = 0;        % fluctuations on hidden states
DCM.options.nonlinear  = 0;        % interactions among hidden states
DCM.options.centre     = 0;        % mean-centre inputs

% priors on connectivity
%--------------------------------------------------------------------------
DCM.a{1}        = [0 0 0;1 0 0;0 1 0];
DCM.a{2}        = [0 1 0;0 0 1;0 0 0];
DCM.b{1}(:,:,1) = [0 0 0;0 0 0;0 0 0];
DCM.b{1}(:,:,2) = [0 0 0;1 0 0;0 0 0];
DCM.b{1}(:,:,3) = [0 0 0;0 1 1;0 0 0];
DCM.b{2}(:,:,1) = [0 0 0;0 0 0;0 0 0];
DCM.b{2}(:,:,2) = [0 0 0;0 0 0;0 0 0];
DCM.b{2}(:,:,3) = [0 0 0;0 1 0;0 0 0];
DCM.c           = [1 0 0;0 0 0;0 0 0];

DCM = spm_dcm_fmri_nmm(DCM);