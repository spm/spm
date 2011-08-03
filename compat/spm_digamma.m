function [y] = spm_digamma(x)
% Digamma function (logarithmic derivative of the gamma function)
% FORMAT [y] = spm_digamma(x)
%
% x - nonnegative, real values
% y - gamma function evaluated at each value x
%
%                    digamma(x) = d(log(gamma(x)))/dx
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_digamma.m 4418 2011-08-03 12:00:13Z guillaume $

persistent runonce
if isempty(runonce)
    warning('spm_digamma is deprecated. Use PSI instead.');
    runonce = 1;
end

y = psi(x);
