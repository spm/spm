function [pE,pC] = spm_ssr_priors(pE,pC)
% augments prior moments of a neural mass model for ssr analyses
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
% pC - prior (co)variances
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
% $Id: spm_ssr_priors.m 4402 2011-07-21 12:37:24Z karl $
 
% catch
%--------------------------------------------------------------------------
try, pE.L; catch, pE.L = 1; end

% number of LFP channels
%--------------------------------------------------------------------------
if size(pE.L,1) == 1, n = size(pE.L,2); else, n = 1; end
if size(pE.C,1),      m = size(pE.C,2); else, m = 1; end

% add prior on spectral density of innovations (pink and white coeficients)
%--------------------------------------------------------------------------
pE.a = sparse(2,m); pC.a = sparse(2,m) + 1/16; % neuronal innovations
pE.b = sparse(2,1); pC.b = sparse(2,1) + 1/16; % channel noise non-specific
pE.c = sparse(2,n); pC.c = sparse(2,n) + 1/16; % channel noise specific

% neuronal innovations
%--------------------------------------------------------------------------
pE.d = sparse(4,m); pC.d = sparse(4,m) + 1/32; 





 
