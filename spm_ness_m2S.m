function [Sp] = spm_ness_m2S(m,C)
% Conditional moments of a Gaussian density (polynomial parameterisation)
% FORMAT [p0,X,F,f,NESS] = spm_ness_hd(M,x)
%--------------------------------------------------------------------------
% m  - (Conditional) mean
% C  - (Conditional) covariance
%
% Sp - Polynomial coefficients or parameters of log density
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


P   = inv(-C);                             % covariance
E   = -P*m(:);                             % mean
n   = numel(m);

% order of expansion (o)
%--------------------------------------------------------------------------
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
Sp    = zeros(size(o,2),1);
for i = 1:n
    k     = (sum(o) == 1) & o(i,:) == 1;
    Sp(k) = E(i);
    for j = i:n
        if j == i
            k = (sum(o) == 2) & (o(i,:) == 2);
            Sp(k) = P(i,j);
        else
            k = (sum(o) == 2) & (o(i,:) == 1) & (o(j,:) == 1);
            Sp(k) = P(i,j);
        end
    end
end
