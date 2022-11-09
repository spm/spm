function Sp = spm_ness_N2Sp(m,C,K)
% Convert a Gaussian density into polynomial potential parameters  
% FORMAT Sp = spm_ness_N2Sp(m,C,[K])
%--------------------------------------------------------------------------
% m  - (Gaussian) mean
% C  - (Gaussian) covariance
% K  - Order of polynomial expansion (K = 3 corresponds to quadratic)

% Sp - Polynomial coefficients or parameters of log density
% n  - Dimensionality of state space
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% defaults
%--------------------------------------------------------------------------
if nargin < 3, K = 3; end

% order of expansion (o)
%--------------------------------------------------------------------------
n     = numel(m);
o     = (1:K) - 1;
for i = 2:n
    o = repmat(o,1,K);
    o = [o; kron((1:K) - 1,ones(1,K^(i - 1)))];
end
k     = sum(o) < K;
o     = o(:,k);

% evaluate expectation and precision
%--------------------------------------------------------------------------
P     = -inv(C);                          % negative precision
E     = -P*m(:);                             % mean
Sp    = zeros(size(o,2),1,'like',m);
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
