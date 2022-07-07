function f = spm_Dpdf(x,a)
% Probability Density Function (PDF) of Dirichlet distribution
% FORMAT f = spm_Dpdf(x,a)
% 
% x - Dirichlet variate
% a - Dirichlet parameters (a>0)
% f - PDF of Dirichlet-distribution at point x
%__________________________________________________________________________
%
% spm_Dpdf implements the Probability Density Function for Dirichlet 
% distribution.
%
% Definition:
%--------------------------------------------------------------------------
% See http://en.wikipedia.org/wiki/Dirichlet_distribution
%
% Algorithm:
%--------------------------------------------------------------------------
% Direct computation using logs and MATLAB's implementation of the log of 
% the gamma function (gammaln).
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


%-Check enough arguments
%--------------------------------------------------------------------------
if nargin<2, error('Insufficient arguments'), end

%-Computation
%--------------------------------------------------------------------------
a = a(:);
x = x(:);

f = exp( gammaln(sum(a)) + sum((a-1).*log(x+eps)) - sum(gammaln(a)) );
