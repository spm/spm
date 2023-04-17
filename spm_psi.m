function [A] = spm_psi(a)
% llog normalised Dirichlet probability matrix (columns)
% FORMAT [A] = spm_psi(a)
%
% a  - Dirichlet parameter tensor
%
% This routine returns log(spm_dir_norm(a)). More formally, it
%  approximates the expectation of the log marginals: E[log(X)]: X(i)
% ~ Beta(a(i),a0 - a(i)), for a > 1. See also: psi.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


% log of an expected probability encoded by Dirichlet parameters
%-------------------------------------------------------------------------- 
A  = spm_log(spm_dir_norm(a));

return

% expectation of a log probability encoded by Dirichlet parameters
%-------------------------------------------------------------------------- 
A  = minus(psi(a),psi(sum(a,1)));
A  = max(A,-32);