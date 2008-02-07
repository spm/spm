function [L] = spm_MH_reml_likelihood(h,Y,M);
% likelihood function for spm_MH_reml
% FORMAT [L] = spm_MH_reml_likelihood(h,Y,M);
%
% h - hyperparameters
% Y - residual covariance
%
% L - likelihood p(Y,P)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MH_reml_likelihood.m 1143 2008-02-07 19:33:33Z spm $

% likelihood function
%--------------------------------------------------------------------------
n = length(Y);
m = length(M.Q);

% compute current estimate of covariance
%--------------------------------------------------------------------------
C     = sparse(n,n);
for i = 1:m
    if M.OPT
        C = C + M.Q{i}*exp(h(i));
    else
        C = C + M.Q{i}*h(i);
    end
end
iC    = inv(C);

% log p(Y,P|M)
%--------------------------------------------------------------------------
e   = h - M.hE;
L   = ...
    - trace(iC*Y)/2 ...
    - e'*M.hP*e/2 ...
    - M.N*n*log(2*pi)/2 ...
    - m*log(2*pi)/2 ...
    - M.N*spm_logdet(C)/2 ...
    + spm_logdet(M.hP)/2;

% p(Y,P|M)
%--------------------------------------------------------------------------
L   = exp(L);

