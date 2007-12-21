function [Q] = spm_Q(a,n,q)
% returns an (n x n) autocorrelation matrix for an AR(p) process
% FORMAT [Q] = spm_Q(a,n,q)
%
% A  - vector pf p AR coeficients
% n  - size of Q
% q  - swtich to return precision [default q = 0]
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
% then trace(cov(y)) ~= n.
% The reason the diagonals of corr(y) are not constant is that we 
% are modeling finite length AR sequences, which incur boundary effects 
% at the beginning and end of the sequence.  These are resolved with a
% Toeplitz form
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_Q.m 1043 2007-12-21 20:34:45Z karl $

% default
%--------------------------------------------------------------------------
try, q; catch, q = 0; end


% compute P
%--------------------------------------------------------------------------
p    = length(a);
if q

    % compute P
    %----------------------------------------------------------------------
    A    = [-a(1) (1 + a(1)^2) -a(1)];
    Q    = spdiags(ones(n,1)*A,[-p:p],n,n);

else

    % compute Q
    %----------------------------------------------------------------------
    A    = [1 -a(:)'];
    P    = spdiags(ones(n,1)*A,-[0:p],n,n);
    K    = inv(P);
    K    = K.*(abs(K) > 1e-4);
    Q    = K*K';
    Q    = toeplitz(Q(:,1));

end


