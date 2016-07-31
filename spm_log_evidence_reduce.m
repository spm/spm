function [F,sE,sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC,TOL)
% Return the log-evidence of a reduced model (under Laplace approximation)
% FORMAT [F,sE,sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC)
%
% qE,qC   - posterior expectation and covariance of full model
% pE,pC   - prior expectation and covariance of full model
% rE,rC   - prior expectation and covariance of reduced model
%
% F       - reduced log-evidence: ln p(y|reduced model) - ln p(y|full model)
% [sE,sC] - posterior expectation and covariance of reduced model
%__________________________________________________________________________
% 
% This routine assumes the reduced model is nested within a full model and
% that the posteriors (and priors) are Gaussian. Nested here means that the
% prior precision of the reduced model, minus the prior precision of the
% full model is positive definite. We additionally assume that the prior
% means are unchanged. The two input argument formats are for use with
% spm_argmax.
%
% This version is the same as spm_log_evidence but performs an
% eigen-reduction of the prior covariance matrix to eliminate fixed
% mixtures of parameters (to ensure well conditioned matrix inversion)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_log_evidence_reduce.m 6849 2016-07-31 12:34:33Z karl $
 

% Compute reduced log-evidence
%==========================================================================
 
% check to see if prior oovaiances are structures
%--------------------------------------------------------------------------
if isstruct(pC), pC = diag(spm_vec(pC)); end
if isstruct(rC), rC = diag(spm_vec(rC)); end
 
% fix tolerance for matrix inversions
%--------------------------------------------------------------------------
if nargin < 7, TOL = 1e-8; end
 
% Remove (a priori) null space
%--------------------------------------------------------------------------
RE    = rE;
SE    = qE;
U     = spm_svd(pC,1e-6);
qE    = U'*spm_vec(qE);
pE    = U'*spm_vec(pE);
rE    = U'*spm_vec(rE);
qC    = U'*qC*U;
pC    = U'*pC*U;
rC    = U'*rC*U;
 
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
    pE = spm_vec(RE);
    rE = sC*sE;
    sE = U*rE + pE - U*U'*pE;
    sC = U*sC*U';
    sE = spm_unvec(sE,SE);
end
