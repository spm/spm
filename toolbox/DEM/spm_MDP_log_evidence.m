function F = spm_MDP_log_evidence(qA,pA,rA)
% Bayesian model reduction for Dirichlet hyperparameters
% FORMAT F = spm_MDP_log_evidence(qA,pA,rA)
%
% qA  - sufficient statistics of posterior of full model
% pA  - sufficient statistics of prior of full model
% rA  - sufficient statistics of prior of reduced model
%
% This routine compute the negative log evidence of a reduced model of a
% categorical distribution parameterised in terms of Dirichlet
% hyperparameters (i.e., concentration parameters encoding probabilities).
% It uses Bayesian model reduction to evaluate the evidence for models with
% and without a particular parameter.
% 
% It is assumed that all the inputs are column vectors.
%
% A demonstration of the implicit pruning can be found at the end of this
% routine
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_log_evidence.m 6756 2016-03-25 09:49:08Z karl $


% change in free energy or log model evidence
%--------------------------------------------------------------------------
dA = qA + rA - pA;
F  = spm_betaln(qA) + spm_betaln(rA) - spm_betaln(pA) - spm_betaln(dA);

return

% notes: synaptic homoeostasis in terms of Bayesian model with production.
% This considers a competition between two inputs:
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

subplot(2,2,1), imagesc(x,x,F + (F > 0)*4), axis square xy,
subplot(2,2,2), plot(x,F'),
xlabel('concentration parameter'), ylabel('log evidence'), axis square




