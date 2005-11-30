function [P,S] = spm_DEM_P(n,r,s,dt,pt)
% returns restricted precisions for DEM
% FORMAT [P,S] = spm_DEM_P(n,r,s,dt,pt)
%__________________________________________________________________________
% n    - embedding order
% r    - restriction order
% s    - temporal smoothness - s.d. of kernel {seconds}
% dt   - time interval {seconds}
% pt   - prediction interval {seconds}
%
% P    - (n x n)     Restricted precision P(pt):     L'*inv(L*R*L')*L;
% S    - {i}(n x n)  Adjusted precisions:        L[-d]'*inv(L*R*L')*L[-i]
%
% R    - (n x n)     covariance of e[:](0):                   R
% L    - (r x n)     prediction operator:                     L
%
%                         e[:] <- E*e(t)
%                  <e[:]*e[:]'> = R = E*V*E'
%                  e[:](t + pt) = L(pt)*e[:](t)
%  <e[:](t + pt)*e[:]'(t + pt)> = L(pt)*R*L(pt)'
%                         F(pt) = -e[:]'*L(pt)'*Q*L(pt)*e[:]/2
%                               = -e[:]'*P(pt)*e[:]/2
%                 dF/dv{d}dv{i} = -de[:]/dv'*S{i}*de[:]/dv
%                          S{i} = L[-d](pt)'*Q*L[-i](pt)
%                            Q  = inv(L(pt)*R*L(pt)')
%
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% default prediction interval
%--------------------------------------------------------------------------
if nargin < 5
    pt = spm_DEM_t(n,r,s,dt);
end

% covariance (R) of error derivatives
%--------------------------------------------------------------------------
R      = spm_DEM_R(n,s,dt);

% precision (P) of errors
%--------------------------------------------------------------------------
j      = 1:(n + n);
L      = triu(toeplitz([1 pt.^j./cumprod(j)]));
Q      = inv(L(1:r,1:n)*R*L(1:r,1:n)');
P      = L(1:r,1:n)'*Q*L(1:r,1:n);
for  j = 1:n
    S{j} = L(1:r,[1:n] + r)'*Q*L(1:r,[1:n] + j - 1);
end

