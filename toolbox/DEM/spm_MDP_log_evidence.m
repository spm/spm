function [F,sA,dFdA] = spm_MDP_log_evidence(qA,pA,rA)
% Bayesian model reduction for Dirichlet hyperparameters
% FORMAT [F,sA,dFdA] = spm_MDP_log_evidence(qA,pA,rA)
%
% qA  - sufficient statistics of posterior of full model
% pA  - sufficient statistics of prior of full model
% rA  - sufficient statistics of prior of reduced model
%
% F    - free energy or (negative) log evidence of reduced model
% sA   - sufficient statistics of reduced posterior
% dFdA - total (negative) free energy gradients with respect to rA
%
% This routine computes the negative log evidence of a reduced model of a
% categorical distribution parameterised in terms of Dirichlet
% hyperparameters (i.e., concentration parameters encoding probabilities).
% It uses Bayesian model reduction to evaluate the evidence for models with
% and without a particular parameter.
%
% A demonstration of the implicit pruning can be found at the end of this
% routine
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% change in free energy or log model evidence
%--------------------------------------------------------------------------
sA = qA + rA - pA;
sA = max(sA,0);
F  = spm_betaln(qA) + spm_betaln(rA) - spm_betaln(pA) - spm_betaln(sA);

if nargout < 3, return, end

% dEdA = d/drA (spm_betaln(rA) - spm_betaln(sA))
%--------------------------------------------------------------------------
% spm_betaln = sum(gammaln(x)) - gammaln(sum(x))
% gammaln    = log(gamma(x))
% psi(x)     = digamma(x) = d(log(gamma(x)))/dx = (d(gamma(x))/dx)/gamma(x).
%--------------------------------------------------------------------------
d_betaln = @(x)psi(x) - psi(sum(x));
dFdA     = d_betaln(rA) - d_betaln(sA);


return

% notes: Illustration of synaptic homoeostasis in terms of Bayesian model
% reduction that considers a competition between two inputs:
%--------------------------------------------------------------------------
x     = linspace(1,32,128);
pA    = [1; 1];
rA    = pA;
rA(2) = 8;
for i = 1:numel(x)
    for j = 1:numel(x)
        qA = [x(i);x(j)];
        F(i,j) = spm_MDP_log_evidence(qA,pA,rA);
    end
end

subplot(2,2,1), imagesc(x,x,F + (F > 0)*4),
title('Free energy landscape','FontSize',16), axis square xy
xlabel('concentration parameter'), ylabel('concentration parameter')
subplot(2,2,2), plot(x,F'),  title('log evidence','FontSize',16)
xlabel('concentration parameter'), ylabel('Log-evidence'), axis square




