function [m,C] = spm_ness_cond(n,K,Sp,ni,x)
% Conditional moments of a Gaussian density (polynomial parameterisation)
% FORMAT [p0,X,F,f,NESS] = spm_ness_hd(M,x)
%--------------------------------------------------------------------------
% n  - Dimensionality of state space
% K  - Order of polynomial expansion (K = 3 corresponds to quadratic)
% Sp - Polynomial coefficients or parameters of log density
% 
% ni - States on which to condition (Optional)
% x  - Values of states [default: 0]
%
% m  - (Conditional) mean
% C  - (Conditional) covariance
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness_cond.m 8072 2021-02-28 16:27:02Z karl $


% order of expansion o
%--------------------------------------------------------------------------
o     = (1:K) - 1;
for i = 2:n
    o = repmat(o,1,K);
    o = [o; kron((1:K) - 1,ones(1,K^(i - 1)))];
end
k     = sum(o) < K;
o     = o(:,k);

% evaluate expectation and precision
%--------------------------------------------------------------------------
E     = zeros(n,1);
P     = zeros(n,n);
for i = 1:n
    k     = (sum(o) == 1) & o(i,:) == 1;
    E(i)  = Sp(k);
    for j = i:n
        if j == i
            k = (sum(o) == 2) & (o(i,:) == 2);
            P(i,j) = 2*Sp(k);
        else
            k = (sum(o) == 2) & (o(i,:) == 1) & (o(j,:) == 1);
            P(i,j) = Sp(k);
            P(j,i) = Sp(k);
        end
        
    end
end
m     = -P\E;
C     = inv(-P);


% Conditional moments if requested
%--------------------------------------------------------------------------
if nargin < 4, return,   end
if nargin < 5, x = ni*0; end

mi     = 1:n;
mi(ni) = [];
m1     = m(mi);
m2     = m(ni);
C11    = C(mi,mi);
C12    = C(mi,ni);
C21    = C(ni,mi);
C22    = C(ni,ni);
m      = m1 + C12*(C22\(x - m2));
C      = C11 - C12*(C22\C21);

return

