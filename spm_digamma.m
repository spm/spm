function y = spm_digamma(x)
% Digamma function (logarithmic derivative of the gamma function)
% 
% FORMAT: y = spm_digamma(x)
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
% %W% Will Penny %E%

error('spm_digamma.c has not been compiled')
