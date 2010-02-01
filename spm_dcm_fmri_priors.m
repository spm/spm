function [pE,pC] = spm_dcm_fmri_priors(A,B,C,varargin)
% Returns the priors for a two-state DCM for fMRI.
% FORMAT:
%    for bi-linear DCM: [pE,pC] = spm_dcm_fmri_priors(A,B,C)
%    for nonlinear DCM: [pE,pC] = spm_dcm_fmri_priors(A,B,C,D)
%    for two-state DCM: [pE,pC] = spm_dcm_fmri_priors(A,B,C,D,'2s')
%
% INPUT:
%    A,B,C,D - constraints on connections (1 - present, 0 - absent)
%
% OUTPUT:
%    pE     - prior expectations (connections and hemodynamic)
%    pC     - prior covariances  (connections and hemodynamic)
%__________________________________________________________________________
%
% References for state equations:
% 1. Marreiros AC, Kiebel SJ, Friston KJ. Dynamic causal modelling for
%    fMRI: a two-state model.
%    Neuroimage. 2008 Jan 1;39(1):269-78.
%
% 2. Stephan KE, Kasper L, Harrison LM, Daunizeau J, den Ouden HE,
%    Breakspear M, Friston KJ. Nonlinear dynamic causal models for fMRI.
%    Neuroimage 42:649-662, 2008.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_priors.m 3705 2010-02-01 20:51:28Z karl $



% number of regions
%--------------------------------------------------------------------------
n     = length(A);       % number of regions

% varargin (D for nonlinear coupling)
%--------------------------------------------------------------------------
if nargin > 3, D = varargin{1}; else, D = zeros(n,n,0); end
if nargin > 4, two_states = 1;  else, two_states = 0;   end

% connectivity priors
%==========================================================================

if two_states
    
    % enforce optimisation of intrinsic (I to E) connections
    %----------------------------------------------------------------------
    A     = (A + eye(n,n)) > 0;

    % prior expectations and variances
    %----------------------------------------------------------------------
    pE.A  =  A*32 - 32;
    pE.B  =  B*0;
    pE.C  =  C*0;
    pE.D  =  D*0;

    % prior covariances
    %----------------------------------------------------------------------
    pC.A  =  A/8;
    pC.B  =  B/8;
    pC.C  =  C*8;
    pC.D  =  D/8;

else

    % prior covariances
    %----------------------------------------------------------------------
    A     =  A - diag(diag(A));
    pC.A  =  A/32;
    pC.B  =  B/32;
    pC.C  =  C;
    pC.D  =  D/32;

    % prior expectations
    %----------------------------------------------------------------------
    pE.A  = -speye(n,n);
    pE.B  =  B*0;
    pE.C  =  C*0;
    pE.D  =  D*0;

end

% and add hemodynamic priors
%==========================================================================
pE.transit = sparse(n,1);
pE.decay   = sparse(n,1);
pE.epsilon = sparse(1,1);

pC.transit = sparse(n,1) + exp(-8);
pC.decay   = sparse(n,1) + exp(-8);
pC.epsilon = sparse(1,1) + exp(-8);

pC         = diag(spm_vec(pC));

return
