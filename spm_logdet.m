function [H] = spm_logdet(C)
% returns the log of the determinant of the positive definite matrix C
% FORMAT [H] = spm_logdet(C)
% H = log(det(C))
%
% spm_logdet is a computationally efficient operator that can deal with
% sparse matrices
%___________________________________________________________________________
% Karl Friston $Id$

% assume diagonal form
%---------------------------------------------------------------------------
H     = sum(log(diag(C)));

% invoke det if non-diagonal
%---------------------------------------------------------------------------
[i j] = find(C);
n     = length(C);
if any(i ~= j)
      a = exp(H/n);
      H = H + log(det(C/a));           
end

% invoke svd is rank deficient
%---------------------------------------------------------------------------
if imag(H) | isinf(H)
    s  = svd(full(C));
    H  = sum(log(s(s > n*max(s)*eps)));
end
