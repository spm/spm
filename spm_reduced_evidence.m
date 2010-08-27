function [F] = spm_reduced_evidence(varargin)
% Return the log-evidence of a reduced model (under Laplace approximation)
% FORMAT [F] = spm_reduced_evidence(qE,qC,pE,pC,rE,rC)
% FORMAT [F] = spm_reduced_evidence(qE,qC,pE,pC,priorfun,varargin)
%
% qE,qC - posterior expectation and covariance of full model
% pE,pC - prior expectation and covariance of full model
% rE,rC - prior expectation and covariance of reduced model
%
% F     - reduced log-evidence: ln p(y|reduced model) - ln p(y|full model)
%
%--------------------------------------------------------------------------
% This routine assumed the reduced model is nested within a full model and
% that the posteriors (and priors) are Gaussian. Nested here means that the
% prior precision of the reduced model, minus the prior precision of the
% full model is positive definite. We additionally assume that the prior
% means are unchanged. The two input argument formats are for use with
% spm_argmin.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_reduced_evidence.m 4053 2010-08-27 19:24:12Z karl $
 
% Compute reduced log-evidence
%==========================================================================

% check to see if priros are specifed by a function
%--------------------------------------------------------------------------
qE = varargin{1};
qC = varargin{2};
pE = varargin{3};
pC = varargin{4};
try
    priors = varargin{5}(varargin{6:end});
    rE     = priors{1};
    rC     = priors{2};
catch
    rE     = varargin{5};
    rC     = varargin{6};
end

% tolerance to handle zero variances
%--------------------------------------------------------------------------
TOL = exp(-32);
 
% preliminaries
%--------------------------------------------------------------------------
qE  = spm_vec(qE);
pE  = spm_vec(pE);
rE  = spm_vec(rE);
qP  = spm_inv(qC,TOL);
pP  = spm_inv(pC,TOL);
rP  = spm_inv(rC,TOL);
dP  = qP + rP - pP;
dE  = qP*qE + rP*rE - pP*pE;
dC  = spm_inv(dP,TOL);
 
% log-evidence
%--------------------------------------------------------------------------
F   = spm_logdet(rP*qP*dC*pC) ...
    - (qE'*qP*qE + rE'*rP*rE - pE'*pP*pE - dE'*dC*dE);
F   = F/2;

