function [pE,pC] = spm_ssr_priors(pE,pC)
% augments prior moments of a neural mass model for ssr analyese
% FORMAT [pE,pC] = spm_ssr_priors(pE,pC)
%
% pE - prior expectation
%
% adds
%
% input and noise parameters
%--------------------------------------------------------------------------
%    pE.a - amplitude of AR  component
%    pE.b - amplitude of IID component
%    pE.c - amplitude of AR  noise (channel specific and non-specific)
%    pE.d - amplitude of IID noise (channel specific and non-specific)
%
%--------------------------------------------------------------------------
%
% pC - prior covariances: cov(spm_vec(pE))
%
% Because priors are specified under log normal assumptions, most
% parameters are simply scaling coefficients with a prior expectation
% and variance of one.  After log transform this renders pE = 0 and
% pC = 1;  The prior expectations of what they scale are specified in
% spm_lfp_fx
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ssr_priors.m 3517 2009-10-29 15:11:56Z guillaume $
 

% add prior on endogenous inputs (neuronal) and noise
%--------------------------------------------------------------------------
pE.a  = 0;               pV.a = 1/16;              % amplitude input AR
pE.b  = 0;               pV.b = 0;                 % amplitude input IID
pE.c  = [0 0];           pV.c = [1/16 1/16];       % amplitude noise AR
pE.d  = [0 0];           pV.d = [1/16 1/16];       % amplitude noise IID

% and augment prior covariance
%--------------------------------------------------------------------------
pC    = spm_cat(spm_diag({pC, diag(spm_vec(pV))}));
 
