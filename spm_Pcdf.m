function F = spm_Pcdf(x,l)
% Cumulative Distribution Function (PDF) of Poisson distribution
% FORMAT F = spm_Ppdf(x,l)
%
% x - ordinates
% l - Poisson mean parameter (lambda l>0) [Defaults to 1]
% F - Poisson CDF
%__________________________________________________________________________
%
% spm_Pcdf implements the Cumulative Distribution Function of the
% Poisson distribution.
% 
% Definition:
%--------------------------------------------------------------------------
% The Poisson Po(l) distribution is the distribution of the number of
% events in unit time for a stationary Poisson process with mean
% parameter lambda=1, or equivalently rate 1/l. If random variable X is
% the number of such events, then X~Po(l), and the CDF F(x) is
% Pr({X<=x}.
%
% F(x) is defined for strictly positive l, given by: (See Evans et al., Ch31)
%
%           { 0                                         for x<0
%           |    _ floor(x)
%    f(rx = |    >      l^i * exp(-l) / i!)             for x>=0
%           {    - i=0
%
% Algorithm:
%--------------------------------------------------------------------------
% F(x), the CDF of the Poisson distribution, for X~Po(l), is related
% to the incomplete gamma function, by:
%
%    F(x) = 1 - gammainc(l,x+1) (x>=0)
%
% See Press et al., Sec6.2 for further details.
%
% Normal approximation:
%--------------------------------------------------------------------------
% For large lambda the normal approximation Y~:~N(l,l) may be used.
% With continuity correction this gives
% F(x) ~=~ Phi((x+.5-l)/sqrt(l))
% where Phi is the standard normal CDF, and ~=~ means "appox. =".
%
% References:
%--------------------------------------------------------------------------
% Evans M, Hastings N, Peacock B (1993)
%       "Statistical Distributions"
%        2nd Ed. Wiley, New York
%
% Abramowitz M, Stegun IA, (1964)
%       "Handbook of Mathematical Functions"
%        US Government Printing Office
%
% Press WH, Teukolsky SA, Vetterling AT, Flannery BP (1992)
%       "Numerical Recipes in C"
%        Cambridge
%
%__________________________________________________________________________
% Copyright (C) 1996-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_Pcdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<2, l=1; end
if nargin<1, error('Insufficient arguments'), end

ad = [ndims(x);ndims(l)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(l),ones(1,rd-ad(2))];];
rs = max(as);
xa = prod(as,2)>1;
if all(xa) && any(diff(as(xa,:)))
    error('non-scalar args must match in size');
end


%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
F = zeros(rs);

%-Only defined for l>0. Return NaN if undefined.
md = ( ones(size(x))  &  l>0 );
if any(~md(:))
    F(~md) = NaN;
    warning('SPM:outOfRangePoisson',...
        'Returning NaN for out of range arguments');
end

%-Non-zero only where defined and x>=0
Q  = find( md  &  x>=0 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Ql=Q; else Ql=1; end

%-Compute
F(Q) = 1 - gammainc(l(Ql),x(Qx)+1);
