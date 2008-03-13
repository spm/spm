function X = spm_pinv(A)
% psudoinverse for sparse matrices
% FORMAT X = spm_pinv(A)
%
% X   - matrix
%
% serial orthogonalisation starting with the first column
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_pinv.m 1207 2008-03-13 20:57:56Z karl $
 
% check A 
%--------------------------------------------------------------------------
[m,n] = size(A);
if isempty(A), X = sparse(n,m); return, end


if m >= n
end


% pseudoinverse
%--------------------------------------------------------------------------
if n > m
   X = spm_pinv(A')';
else
   [U,S,V] = spm_svd(A,0);
   if m > 1, s = full(diag(S));
      elseif m == 1, s = S(1);
      else s = 0;
   end
   tol = max(m,n)*eps(max(s));
   r   = sum(s > tol);
   if ~r
      X = sparse(n,m);
   else
      i = 1:r;
      s = sparse(i,i,1./s(i),r,r);
      X = V(:,i)*s*U(:,i)';
   end
end
