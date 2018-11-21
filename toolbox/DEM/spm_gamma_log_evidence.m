function [F,sA] = spm_gamma_log_evidence(qA,pA,rA)
% Bayesian model reduction for gamma distibutions
% FORMAT [F,sA] = spm_gamma_log_evidence(qA,pA,rA)
%
% qA  - 2-vector with shape/rate parameter of posterior of full model
% pA  - 2-vector with shape/rate parameter of prior of full model
% rA  - 2-vector with shape/rate parameter of prior of reduced model
%
%
% F   - (negative) free energy or log evidence of reduced model
% sA  - shape/rate parameter of reduced posterior
%
% This routine computes the negative log evidence of a reduced model of a
% gamma distribution parameterised in terms of its shape parameter.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gamma_log_evidence.m 7493 2018-11-21 12:47:55Z thomas $

% reduced posteriors
%--------------------------------------------------------------------------
sA = qA + rA - pA;

% change in free energy or log model evidence - shape parameter
%--------------------------------------------------------------------------
F  = - gammaln(qA(1)) - gammaln(rA(1)) + gammaln(pA(1)) + gammaln(sA(1));

% change in free energy or log model evidence - rate parameter
%--------------------------------------------------------------------------
F  = F + qA(1)*log(qA(2)) + rA(1)*log(rA(2)) - pA(1)*log(pA(2)) - sA(1)*log(sA(2));
