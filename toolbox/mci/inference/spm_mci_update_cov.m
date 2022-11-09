function [P] = spm_mci_update_cov (P)
% Update covariance matrix of proposal density using Robbins-Monro
% FORMAT [P] = spm_mci_update_cov (P)
%
% See e.g.
% H. Haario, E. Saksman, and J. Tamminen. An adaptive Metropolis algorithm. 
% Bernoulli, 7(2):223-242, 2001.
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

Np=size(P.theta,1);
gamma=1/P.adapt_its;

x=P.theta(:,end);
dx=x-P.mu;

P.mu=P.mu+gamma*dx;
P.Ct=P.Ct+gamma*(dx*dx'-P.Ct);

