function [A] = spm_psi(A)
% Normalisation of a probability transition rate matrix (columns)
% FORMAT [A] = spm_psi(A)
%
% A  - numeric array
%
% See also: psi.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


% normalisation of a probability transition rate matrix (columns)
%--------------------------------------------------------------------------
A = bsxfun(@minus, psi(A), psi(sum(A,1)));
