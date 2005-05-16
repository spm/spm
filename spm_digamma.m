function [y] = spm_digamma(x)
% Digamma function (logarithmic derivative of the gamma function)
% FORMAT [y] = spm_digamma(x)
%
% x - nonnegative, real values
% y - gamma function evaluated at each value x
%
%                    digamma(x) = d(log(gamma(x)))/dx
%
% Reference:
%   Mex file derived from a FORTRAN program by D. E. Amos,
%   ACM Transactions on Mathematical Software, 1983.
%   Obtained from NETLIB: http://www.netlib.org/toms/610.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny
% $Id: spm_digamma.m 159 2005-05-16 14:00:56Z guillaume $

error('spm_digamma.c not compiled - see Makefile')
