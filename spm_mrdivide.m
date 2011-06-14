function D = spm_mrdivide(A, B)
% Regularised variant of mrdivide(A, B) or A / B, similar to B * spm_inv(A)
% FORMAT D = spm_mrdivide(A, B)
%
% D = B * inv(A), or if A is near singular D = B * inv(A + TOL*eye(size(A))
% 
% where TOL is adaptively increased if necessary.
%
% This function should be preferable to B * spm_inv(A) if A is large and
% sparse or if B has few rows, since the inverse need not be explicitly
% computed (the linear system can be solved with the backslash operator).
%
% See also: spm_mldivide
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Ged Ridgway
% $Id: spm_mrdivide.m 4360 2011-06-14 16:46:37Z ged $

D = spm_mldivide(B', A')';
