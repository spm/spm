function A = spm_log(A)
% Log of numeric array plus a small constant
% FORMAT A = spm_log(A)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

A          = log(A);
A(A < -32) = -32;