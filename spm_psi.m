function [A] = spm_psi(A)
% normalisation of a probability transition rate matrix (columns)
% FORMAT [A] = spm_psi(A)
%
% A  - numeric array
%
% See also: psi.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_psi.m 7300 2018-04-25 21:14:07Z karl $

% normalisation of a probability transition rate matrix (columns)
%--------------------------------------------------------------------------
A = bsxfun(@minus, psi(A), psi(sum(A,1)));
