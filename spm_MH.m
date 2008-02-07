function [P,F] = spm_MH(L,B,y,M)
% The Rejection-Metropolis-Hastings Algorithm
% FORMAT [P,F] = spm_MH(L,P,y)
%
% L   - likelihood function: inline(P,y,M)
% B   - free parameter [structure]
% Y   - response  [stucture]
% M   - model [structure]
%
% P   - Sample from posterior p(P|y,M)
% F   - marginal likelihood p(y|M) using harmonic mean
%--------------------------------------------------------------------------
%
% Returns a harmonic mean estimate of the log-marginal likelihood or
% log-evidence and a sample from the posterior density of the free parameters
% of a model.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MH.m 1143 2008-02-07 19:33:33Z spm $

% initialise parameters
%--------------------------------------------------------------------------
P(:,1) = spm_vec(B);

% MCMC - RMH
%--------------------------------------------------------------------------
n     = 2^8;                                          % number of burn in
N     = 2^16;                                         % number of samples
for i = 1:N

    % sample from proposal
    %----------------------------------------------------------------------
    pi = P(:,i);
    pp = pi + randn(size(P,1),1)/32;

    % compute importnace ratio
    %----------------------------------------------------------------------
    Lp = feval(L,spm_unvec(pp,B),y,M);
    Li = feval(L,spm_unvec(pi,B),y,M);
    r  = Lp/Li;

    % accept and marginal likelihood
    %----------------------------------------------------------------------
    if rand < r
        P(:,i + 1) = pp;
        F(i + 1)   = Lp;
    else
        P(:,i + 1) = pi;
        F(i + 1)   = Li;
    end

end

% burn in
%--------------------------------------------------------------------------
P = P(:,n:end);
F = F(n:end);
F = -log(mean(mean(F)./F)) + log(mean(F));
