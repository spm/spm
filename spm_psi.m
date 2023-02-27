function [A] = spm_psi(a)
% Normalisation of a Dirichlet probability matrix (columns)
% FORMAT [A] = spm_psi(a)
%
% a  - Dirichlet parameter tensor
%
% This can be regarded as log(spm_dir_norm(a)). More formally, it
% corresponds to  the expectation  of the log marginals: E[log(X)]: X(i)
% ~ Beta(a(i),a0 - a(i)). See also: psi.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


% expectation of a log probability encoded by Dirichlet parameters
%-------------------------------------------------------------------------- 
A  = minus(psi(a),psi(sum(a,1)));
A  = max(A,-32);