function [L,L2,st] = spm_mci_joint (Pr,M,U,Y,beta)
% Compute log joint probability of model
% FORMAT [L,L2,st] = spm_mci_joint (Pr,M,U,Y,beta)
%
% Pr    parameters (vectorised and in M.V subspace)
% M     model structure
% U     inputs
% Y     data
% beta  inverse temperature
%
% L     beta * log p(Y|P) + log p(P)
% L2    log p(Y|P)
% st    status flag (0 for OK, -1 for problem)
%
% A default beta=1 gives usual log joint
%__________________________________________________________________________

% Will Penny and Biswa Sengupta
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

st=0;
if nargin < 5 | isempty(beta)
    beta=1;
end
Pr=Pr(:);

% Parameter errors in subspace
e = Pr;
L1 = - e'*M.ipC*e/2 + M.log_prior_t2;

% Parameters in original space
P = M.V*Pr+M.vpE;

[L2,tmp,st] = feval(M.L,P,M,U,Y);
L = L1+beta*L2;