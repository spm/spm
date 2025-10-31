function [L] = spm_MDP_BMR(qp,rp)
% returns the likelihhod of reduced (Dirichlet) models
% FORMAT [L] = spm_MDP_BMR(qp,rp)
%
% qp    - posterior Dirichlet parameters (tensor or numeric array)
% rp{m} - cell array of m reduced priors
%
% L(m)  - likelihood of model m
%
% This auxillary routine assumes the priors are the (model) average of
% reduced priors.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

% get posteriors and priors
%--------------------------------------------------------------------------
np    = numel(rp);
pp    = 0;
for i = 1:np
    pp = pp + rp{i};
end
pp    = pp/sum(pp,'all');

% Bayesian model reduction
%------------------------------------------------------
F     = zeros(np,1);
for i = 1:np
    F(i) = sum(spm_MDP_log_evidence(qp,pp,rp{i}),'all');
end
L     = spm_softmax(-F);
