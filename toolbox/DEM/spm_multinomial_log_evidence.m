function [F,sA] = spm_multinomial_log_evidence(qA,pA,rA)
% Bayesian model reduction for multinomial distibutions
% FORMAT [F,sA] = spm_multinomial_log_evidence(qA,pA,rA)
%
% qA  - parameter of posterior of full model
% pA  - parameter of prior of full model
% rA  - parameter of prior of reduced model
%
%
% F   - (negative) free energy or log evidence of reduced model
% sA  - parameter of reduced posterior
%
% This routine computes the negative log evidence of a reduced model of a
% mutinomial distribution. This also applies for Bernoulli, Binomial, and
% Categorical distributions.
%__________________________________________________________________________

% Thomas Parr
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% reduced posteriors
%--------------------------------------------------------------------------
sA = spm_softmax(log(qA) + log(rA) - log(pA));

% change in free energy or log model evidence
%--------------------------------------------------------------------------
k = find(rA);
F  = log(rA(k(1))) - log(pA(k(1))) + log(qA(k(1))) - log(sA(k(1)));




