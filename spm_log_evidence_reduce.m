function [F,sE,sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC)
% Return the log-evidence of a reduced model (under Laplace approximation)
% FORMAT [F,sE,sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC)
%
% qE,qC    - posterior expectation and covariance of full model
% pE,pC    - prior expectation and covariance of full model
% rE,rC    - prior expectation and covariance of reduced model
%
% F        - reduced log-evidence: ln p(y|reduced model) - ln p(y|full model)
% [sE,sC]  - posterior expectation and covariance of reduced model
%
%--------------------------------------------------------------------------
% This routine assumes the reduced model is nested within a full model and
% that the posteriors (and priors) are Gaussian. Nested here means that the
% prior precision of the reduced model, minus the prior precision of the
% full model is positive definite. We additionally assume that the prior
% means are unchanged. The two input argument formats are for use with
% spm_argmax.
%
% This version is the same as spm_log_evidence but perform an
% eigen-reduction of the prior covariance matrix to eliminate fixed mixtures
% of parameters ()and ensure well conditioned matrix inversion)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_log_evidence_reduce.m 5022 2012-10-30 19:25:02Z karl $
 
% Compute reduced log-evidence
%==========================================================================
 
% check to see if priors are specified by a function
%--------------------------------------------------------------------------
 
 
% Remove (a priori) null space
%--------------------------------------------------------------------------
E     = pE;
U     = spm_svd(pC);
qE    = U'*spm_vec(qE);
pE    = U'*spm_vec(pE);
rE    = U'*spm_vec(rE);
qC    = U'*qC*U;
pC    = U'*pC*U;
rC    = U'*rC*U;
 
% fix tolerance for matrix inversions
%--------------------------------------------------------------------------
TOL   = exp(-16);
 
% preliminaries
%--------------------------------------------------------------------------
qP    = spm_inv(qC,TOL);
pP    = spm_inv(pC,TOL);
rP    = spm_inv(rC,TOL);
sP    = qP + rP - pP;
sC    = spm_inv(sP,TOL);
pC    = spm_inv(pP,TOL);
sE    = qP*qE + rP*rE - pP*pE;
 
% log-evidence
%--------------------------------------------------------------------------
F     = spm_logdet(rP*qP*sC*pC) ...
      - (qE'*qP*qE + rE'*rP*rE - pE'*pP*pE - sE'*sC*sE);
F     = F/2;
    
% restore full conditional density
%--------------------------------------------------------------------------
if nargout > 1
    rE = sC*sE;
    rC = sC;
    sE = spm_unvec(U*rE,E);
    sC = U*rC*U';
end
