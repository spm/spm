function A = spm_log(A)
% log of numeric array plus a small constant
% FORMAT A = spm_log(A)
%--------------------------------------------------------------------------
A           = log(A);
A(isinf(A)) = -32;