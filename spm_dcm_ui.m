function varargout = spm_dcm_ui(Action)
% User interface for Dynamic Causal Modelling (DCM)
% FORMAT spm_dcm_ui('specify')
% FORMAT spm_dcm_ui('estimate')
% FORMAT spm_dcm_ui('search')
% FORMAT spm_dcm_ui('optimise')
% FORMAT spm_dcm_ui('review')
% FORMAT spm_dcm_ui('compare')
% FORMAT spm_dcm_ui('average (BPA)')
% FORMAT spm_dcm_ui('average (BMA)')
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
%   DCM.Pp     - posterior probabilities
%   DCM.Vp     - variance of parameter estimates
%   DCM.H1     - 1st order Volterra Kernels - hemodynamic
%   DCM.K1     - 1st order Volterra Kernels - neuronal
%   DCM.R      - residuals
%   DCM.y      - predicted responses
%   DCM.xY     - original response variable structures
%   DCM.T      - threshold for inference based on posterior p.d.f
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
% Copyright (C) 2002-2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_ui.m 6005 2014-05-21 16:46:26Z guillaume $


DCMversion = 'DCM12';

if nargin == 1 && strcmpi(Action,'version')
    varargout = {DCMversion};
    return
end

% Get figure handles
%--------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
header = get(Finter,'Name');
spm_clf(Finter);
set(Finter,'Name',sprintf('Dynamic Causal Modelling (%s)',DCMversion));
spm('Pointer','Arrow');

% Temporary welcome message
%--------------------------------------------------------------------------
disp(['Please refer to this version as ' ...
    DCMversion ' in papers and publications.']);

% Options, using pull-down menu
%--------------------------------------------------------------------------
if ~nargin
    str       =  'Action: ';
    Actions   = {'specify',  ...
                 'estimate (time-series)', ...
                 'estimate (cross-spectra)', ...
                 'search', ...
                 'optimise', ...
                 'review',   ...
                 'compare',  ...
                 'average',  ...
                 'quit'};
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
% Estimate models - standard
%==========================================================================
case 'estimate (time-series)',
    
    %-estimate models
    %----------------------------------------------------------------------
    spm('FnBanner','spm_dcm_estimate');

    %-select DCM models
    %----------------------------------------------------------------------
    [P, sts] = spm_select([1 Inf],'^DCM.*\.mat$','select DCM_???.mat');
    if ~sts, return; else P = cellstr(P); end

    spm('Pointer','Watch');
    spm('FigName','Estimation in progress');

    %-loop over models
    %----------------------------------------------------------------------
    for i=1:numel(P)
        spm('SFnBanner',sprintf('spm_dcm_estimate: model %d',i));
        spm_dcm_estimate(P{i});
    end
    
%==========================================================================  
% Estimate models - cross spectral density
%==========================================================================
case 'estimate (cross-spectra)',
    
    %-estimate models
    %----------------------------------------------------------------------
    spm('FnBanner','spm_dcm_fmri_csd');

    %-select DCM models
    %----------------------------------------------------------------------
    [P, sts] = spm_select([1 Inf],'^DCM.*\.mat$','select DCM_???.mat');
    if ~sts, return; else P = cellstr(P); end

    spm('Pointer','Watch');
    spm('FigName','Estimation in progress');

    %-loop over models
    %----------------------------------------------------------------------
    for i=1:numel(P)
        spm('SFnBanner',sprintf('spm_dcm_fmri_csd: model %d',i));
        spm_dcm_fmri_csd(P{i});
    end

%==========================================================================
% Esimate and search a model set
%==========================================================================
case 'search',

    spm('FnBanner','spm_dcm_search');
    
    spm_dcm_search;

    
%==========================================================================
% Post hoc model optimisation/selection
%==========================================================================
case 'optimise',

    spm('FnBanner','spm_dcm_post_hoc');
    
    spm_dcm_post_hoc;

    
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
    
    spm('FnBanner','spm_api_bmc');
    
    spm_jobman('Interactive','','spm.dcm.bms.inference');


%==========================================================================
% Average
%==========================================================================
case 'average',
    
    if spm_input('Average',1,'b',{'BPA','BMA'},[1 0])
        
        spm('FnBanner','spm_dcm_average');
        spm_dcm_average;         %  Average several models (Bayesian FFX)
        
    else
        
        spm('FnBanner','spm_dcm_bma_results');
        spm_dcm_bma_results;     %  Average model parameters from BMS (BMA)
        
    end
    
    
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
spm_input('Thank you',1,'d');
