function [varargout] = spm_dcm_optimise(qE,qC,pE,pC,priorfun,varargin)
% Optimises the priors of a model (under Laplace approximation)
% FORMAT [rE,rC] = spm_dcm_optimise(qE,qC,pE,pC)
%
% qE,qC    - posterior expectation and covariance of model
% pE,pC    - prior expectation and covariance of model
% priorfun - inline function that returns priors
%            {rE rC} = priorfun(varargin{:})
%
% pE,pC    - optimal priors defining a reduced model
%
%--------------------------------------------------------------------------
% This routine assumed the reduced model is nested within a full model and
% that the posteriors (and priors) are Gaussian. Nested here means that the
% prior precision of the reduced model, minus the prior precision of the
% full model is positive definite. We additionally assume that the prior
% means are unchanged. 
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_optimise.m 4053 2010-08-27 19:24:12Z karl $
 
% Compute reduced log-evidence
%==========================================================================

% default prior function
%--------------------------------------------------------------------------
if nargin == 4
    priorfun    = inline('{v2, diag(exp(v1))}','v1','v2');
    varargin{1} = log(diag(pC));
    varargin{2} = pE;
end

varargout   = spm_argmax('spm_reduced_evidence',qE,qC,pE,pC,priorfun,varargin{:},6);
varargin{1} = varargout;
varargout   = priorfun(varargin{:});

