function [Q] = spm_Q(A,n)
% return autocorrrelation mtrix of an n x n AR(p) process
% FORMAT [Q] = spm_Q(A,n)
%
% A  - vector pf p AR coeficients
% n  - size of Q
%___________________________________________________________________________
% %W% Karl Friston %E%

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
