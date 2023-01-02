function [A] = spm_psi(a)
% Normalisation of a Dirichlet probability matrix (columns)
% FORMAT [A] = spm_psi(a)
%
% a  - Dirichlet tensor
%
% This can be regarded as log(spm_dir_norm(a)). More formally, it
% corresponds to  the expectation  of the log marginals: E[log(X)]: X(i)
% ~ Beta(a(i),a0 - a(i)). See also: psi.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


% normalisation of a probability transition rate matrix (columns)
%--------------------------------------------------------------------------
A = minus(psi(a),psi(sum(a,1)));
A(A < -32) = -32;