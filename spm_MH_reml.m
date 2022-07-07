function [F,P] = spm_MH_reml(YY,X,Q,N,hE);
% Estimation of covariance components from y*y' using sampling
% FORMAT [F,P] = spm_MH_reml(YY,X,Q,N,[hE]);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
% N   - number of samples
%
% hE  - prior expectation: log-normal hyper-parameterisation (with hyperpriors)
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q)
% P   - sample of hyperparameters from their posterior p(h|YY,X,Q)
%--------------------------------------------------------------------------
%
% This routine is using MCMC sampling (reverible Metropolis-Hastings)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2006-2022 Wellcome Centre for Human Neuroimaging


% assume a single sample if not specified
%--------------------------------------------------------------------------
try
    N;
catch
    N  = 1;
end

% assume OPT = 0
%--------------------------------------------------------------------------
try
    hE;
    OPT = hE;
catch
    OPT = 0;
end

% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(length(Q{1}),1);
else
    X = orth(full(X));
end

% remove fixed effects
%--------------------------------------------------------------------------
n     = length(Q{1});
m     = length(Q);
h     = zeros(m,1);
R     = speye(n,n) - X*X';
YY    = R*YY*R;

M.OPT = OPT;
M.Q   = Q;
M.N   = N;

% initialise and specify hyperpriors
%--------------------------------------------------------------------------
[C,h,Ph,Fr] = spm_reml(YY,X,Q,N,0,4,OPT);
if M.OPT
    M.hE  = h - 16;
    M.hP  = eye(m,m)/32;
else
    M.hE  = zeros(m,1);
    M.hP  = speye(m,m)/exp(32);
end

% sample
%--------------------------------------------------------------------------
[P,F] = spm_MH('spm_MH_reml_likelihood',h,YY,M);
