function [F,sA] = spm_gamma_log_evidence(qA,pA,rA)
% Bayesian model reduction for Dirichlet hyperparameters
% FORMAT [F,sA] = spm_gamma_log_evidence(qA,pA,rA)
%
% qA  - shape parameter of posterior of full model
% pA  - shape parameter of prior of full model
% rA  - shape parameter of prior of reduced model
%
% F   - (negative) free energy or log evidence of reduced model
% sA  - shape parameter of reduced posterior
%
% This routine compute the negative log evidence of a reduced model of a
% gamma distribution parameterised in terms of its shape parameter.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gamma_log_evidence.m 7326 2018-06-06 12:16:40Z karl $


% change in free energy or log model evidence
%--------------------------------------------------------------------------
sA = qA + rA - pA;
F  = - gammaln(qA) - gammaln(rA) + gammaln(pA) + gammaln(sA);




