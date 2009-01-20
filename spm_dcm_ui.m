function spm_dcm_ui(Action)
% User interface for Dynamic Causal Modelling (DCM)
% FORMAT spm_dcm_ui('specify')
% FORMAT spm_dcm_ui('estimate')
% FORMAT spm_dcm_ui('review')
% FORMAT spm_dcm_ui('compare')
% FORMAT spm_dcm_ui('average')
%
% * Specify a new model
% * Estimate a specified model
% * Review a previously estimated model
% * Compare two or more estimated models
% * Produce an aggregate model using Bayesian averaging
%
% DCM structure, as saved in DCM_???.mat:
%
%   DCM.M      - model  specification structure (see spm_nlsi)
%   DCM.Y      - output specification structure (see spm_nlsi)
%   DCM.U      - input  specification structure (see spm_nlsi)
%   DCM.Ep     - posterior expectations (see spm_nlsi)
%   DCM.Cp     - posterior covariances (see spm_nlsi)
%   DCM.a      - intrinsic connection matrix
%   DCM.b      - input-dependent connection matrix
%   DCM.c      - input connection matrix
%   DCM.pA     - pA - posterior probabilities
%   DCM.pB     - pB - posterior probabilities
%   DCM.pC     - pC - posterior probabilities
%   DCM.vA     - vA - variance of parameter estimates
%   DCM.vB     - vB - variance of parameter estimates
%   DCM.vC     - vC - variance of parameter estimates
%   DCM.H1     - 1st order Volterra Kernels - hemodynamic
%   DCM.H2     - 1st order Volterra Kernels - hemodynamic
%   DCM.K1     - 1st order Volterra Kernels - neuronal
%   DCM.K1     - 1st order Volterra Kernels - neuronal
%   DCM.R      - residuals
%   DCM.y      - predicted responses
%   DCM.xY     - original response variable structures
%   DCM.T      - threshold for inference based on posterior p.d.f
%   DCM.Ce     - Estimated observation noise covariance
%   DCM.v      - Number of scans
%   DCM.n      - Number of regions
%
%__________________________________________________________________________
%
% DCM  is a  causal  modelling  procedure  for dynamical  systems  in which
% causality  is inherent  in the  differential equations  that  specify the
% model.  The basic idea  is to treat the system of interest,  in this case
% the brain,  as an  input-state-output  system.  By perturbing  the system
% with  known  inputs,  measured  responses  are used to  estimate  various
% parameters that govern the evolution of brain states.  Although there are
% no  restrictions  on  the  parameterisation  of  the  model,  a  bilinear
% approximation affords a simple  re-parameterisation in terms of effective
% connectivity.  This effective connectivity can be latent or intrinsic or,
% through  bilinear  terms,  model  input-dependent  changes  in  effective
% connectivity.   Parameter  estimation   proceeds  using  fairly  standard
% approaches to system identification that rest upon Bayesian inference.
% 
% Dynamic  causal  modelling   represents  a   fundamental  departure  from
% conventional   approaches   to   modelling   effective   connectivity  in
% neuroscience.  The critical distinction between DCM and other approaches,
% such as  structural equation  modelling  or  multivariate  autoregressive
% techniques is that the input is treated as known, as opposed to stochastic.
% In this sense DCM is much closer to conventional analyses of neuroimaging
% time series  because the causal  or explanatory  variables enter as known
% fixed quantities.  The use of designed and known inputs in characterising
% neuroimaging data  with the general linear model or DCM is a more natural
% way to  analyse data  from  designed  experiments.  Given  that  the vast
% majority  of imaging  neuroscience  relies  upon designed  experiments we
% consider DCM a potentially useful complement to existing techniques.  
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_ui.m 2630 2009-01-20 17:11:57Z cc $


% Get figure handles
%--------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
header = get(Finter,'Name');
spm_clf(Finter);
set(Finter,'Name','Dynamic Causal Modelling');
spm('Pointer','Arrow');


% Options, using pull-down menu
%--------------------------------------------------------------------------
if ~nargin
    str      = 'Action: ';
    Actions  = {'specify','estimate','review','compare','average','quit'};
    selected = spm_input(str,1,'m',Actions);
    Action   = Actions{selected};
end


switch lower(Action)

%==========================================================================
% Specify graph
%==========================================================================
case 'specify',

    spm('FnBanner','spm_dcm_specify');
    
    spm_dcm_specify;
    
    
%==========================================================================  
% Estimate models
%==========================================================================
case 'estimate',
    
    %-estimate models
    %----------------------------------------------------------------------
    spm('FnBanner','spm_dcm_estimate');

    %-select DCM models
    %----------------------------------------------------------------------
    P = cellstr(spm_select(Inf,'^DCM.*\.mat$','select DCM_???.mat'));

    spm('Pointer','Watch');
    spm('FigName','Estimation in progress');

    %-loop over models
    %----------------------------------------------------------------------
    for i=1:numel(P)
        spm('SFnBanner',sprintf('spm_dcm_estimate: model %d',i));
        spm_dcm_estimate(P{i});
    end

    
%==========================================================================
% Review results
%==========================================================================
case 'review',

    spm('FnBanner','spm_dcm_review');
    
    spm_dcm_review;

    
%==========================================================================
% Compare different models
%==========================================================================
case 'compare',
    
    %spm('FnBanner','spm_dcm_compare');
    
    spm_jobman('Interactive','','spm.stats.bms.bms_dcm')

    %spm_dcm_compare;

    
%==========================================================================
% Average several models (Bayesian FFX)
%==========================================================================
case 'average',
    
    spm('FnBanner','spm_dcm_average');
    
    spm_dcm_average(1);  % ERP: 0; fMRI: any value > 0

    
%==========================================================================
% Quit DCM GUI
%==========================================================================
case 'quit',
    

%==========================================================================
% Otherwise
%==========================================================================
otherwise
    error('Unknown action string.');
    
end

% Return to SPM
%--------------------------------------------------------------------------
spm_clf(Finter);
spm('Pointer','Arrow');
set(Finter,'Name',header);
if strcmpi(Action, 'compare')
    spm_input('Thank you',3,'d');
else
    spm_input('Thank you',1,'d');
end

return
