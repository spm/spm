function [Sp,o] = spm_ness_cond_inv(m,C)
% Returns polynomial parameterisation (sP)of a Gaussian density
% FORMAT [Sp,o] = spm_ness_cond_inv(m,C)
%--------------------------------------------------------------------------
% m  - mean
% C  - covariance
%
% Sp - Polynomial coefficients; i.e., parameters of log density
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% order of expansion (o)
%--------------------------------------------------------------------------
n     = numel(m);
K     = 3;
o     = (1:K) - 1;
for i = 2:n
    o = repmat(o,1,K);
    o = [o; kron((1:K) - 1,ones(1,K^(i - 1)))];
end
k     = sum(o) < K;
o     = o(:,k);


% evaluate expectation and precision
%--------------------------------------------------------------------------
P     = -spm_inv(C);                         % precision
E     = -P*m ;                           % expectation
for i = 1:n
    k     = (sum(o) == 1) & o(i,:) == 1;
    Sp(k,1) = E(i);
    for j = i:n
        if j == i
            k = (sum(o) == 2) & (o(i,:) == 2);
            Sp(k,1) = P(i,j);
        else
            k = (sum(o) == 2) & (o(i,:) == 1) & (o(j,:) == 1);
            Sp(k,1) = P(i,j);
            Sp(k,1) = P(j,i);
        end
    end
end

