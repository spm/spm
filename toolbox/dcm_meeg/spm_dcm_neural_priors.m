function [pE,pC] = spm_dcm_neural_priors(A,B,C,model)
% Prepares the priros on the paramters of neural mass models
% FORMAT [pE,pC] = spm_dcm_neural_priors(A,B,C,'model'))
%
% A,B{m},C  - binary constraints on extrinsic connections for m conditions
% 'model'   - 'ERP','SEP','LFP','NNM' or 'MFM'
%
% pE - prior expectation - f(x,u,P,M)
%
% synaptic parameters
%--------------------------------------------------------------------------
%    pE.T - syaptic time constants
%    pE.H - syaptic densities
%    pE.S - activation function parameters
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A  - extrinsic
%    pE.B  - trial-dependent
%    pE.C  - stimulus input
%
%    pE.SA - switches on extrinsic (excitatory)
%    pE.GE - switches on intrinsic (excitatory)
%    pE.GI - switches on intrinsic (inhibitory)
%
%    pE.D  - delays
%
% stimulus and noise parameters
%--------------------------------------------------------------------------
%    pE.R - onset and dispersion
%    pE.X - endogenous activity
%
% pC - prior covariances: cov(spm_vec(pE))
%
% Because priors are specified under log normal assumptions, most
% parameters are simply scaling coefficients with a prior expectation
% and variance of one.  After log transform this renders pE = 0 and
% pC = 1;
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_neural_priors.m 2374 2008-10-21 18:52:29Z karl $

% check options
%==========================================================================

% get priors on neural model
%--------------------------------------------------------------------------
switch lower(model)

    % linear David et al model (linear in states)
    %======================================================================
    case{'erp'}

        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_erp_priors(A,B,C);


    % linear David et al model (linear in states) - fast version for SEPs
    %======================================================================
    case{'sep'}

        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC]   = spm_erp_priors(A,B,C);
        pE.T      = pE.T - 1;
        pE.R(:,2) = pE.R(:,2) - 2;
        
        
    % linear David et al model (linear in states)  - with self-inhibition
    %======================================================================
    case{'lfp'}

        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_lfp_priors(A,B,C);


    % Neural mass model (nonlinear in states)
    %======================================================================
    case{'nmm','mfm'}

        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_nmm_priors(A,B,C);

    otherwise
        warndlg('Unknown model')
end
