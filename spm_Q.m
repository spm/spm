function [Q] = spm_Q(A,n)
% returns an (n x n) autocorrelation matrix for an AR(p) process
% FORMAT [Q] = spm_Q(A,n)
%
% A  - vector pf p AR coeficients
% n  - size of Q
%___________________________________________________________________________
% @(#)spm_Q.m	2.2 Karl Friston 03/03/03

% compute Q
%---------------------------------------------------------------------------
p    = length(A);
A    = [1 -A(:)'];
K    = inv(spdiags(ones(n,1)*A,-[0:p],n,n));
K    = K.*(abs(K) > 1e-4);
Q    = K*K';
D    = spdiags(sqrt(1./diag(Q)),0,n,n);
Q    = D*Q*D;
Q    = Q.*(abs(Q) > 1e-4);
