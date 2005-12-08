function [Q] = spm_Q(A,n)
% returns an (n x n) autocorrelation matrix for an AR(p) process
% FORMAT [Q] = spm_Q(A,n)
%
% A  - vector pf p AR coeficients
% n  - size of Q
%__________________________________________________________________________
% spm_Q uses a Yule-Walker device to compute K where:
% 
% y = K*z
% 
% such that y is an AR(n) process generated from an i.i.d innovation 
% z.  This means
% 
% cov(y) = <K*z*z'*K> = K*K'
% 
% Critically, this is not the correlation because if cov(z) = eye(n) 
% then trace(cov(y)) ~= n.  This is why the normalization is required
% 
% corr(y) = D*K*K'*D';
% 
% The reason the diagonals of corr(y)  are not constant is that we 
% are modeling finite length AR sequences, which incur boundary effects 
% at the beginning and end of the sequence.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_Q.m 372 2005-12-08 17:12:13Z karl $


% compute Q
%--------------------------------------------------------------------------
p    = length(A);
A    = [1 -A(:)'];
K    = inv(spdiags(ones(n,1)*A,-[0:p],n,n));
K    = K.*(abs(K) > 1e-4);
Q    = K*K';
D    = spdiags(sqrt(1./diag(Q)),0,n,n);
Q    = D*Q*D;
Q    = Q.*(abs(Q) > 1e-4);