function f = spm_Ipdf(x,n,p)
% Probability Distribution Function of Binomial distribution
% FORMAT f = spm_Ipdf(x,n,p)
%
% x - ordinates
% n - Binomial n
% p - Binomial p [Defaults to 0.5]
% f - PDF
%__________________________________________________________________________
%
% spm_Ipdf returns the Probability (Distribution) Function (PDF) for
% the Binomial family of distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The Bin(n,p) distribution is the distribution of the number of
% successes from n identical independent Bernoulli trials each with
% success probability p. If random variable X is the number of
% successes from such a set of Bernoulli trials, then the PDF f(x) is
% Pr{X=x}, defined for whole n (i.e. non-negative integer n) and
% p in [0,1], given by: (See Evans et al., Ch6)
%
%           { nCx * p^x * (1-p)^(n-x)     for x=0,1,...,n
%    f(x) = |
%           {  0                          otherwise
%
% where nCx is the Binomial coefficient "n-choose-x", given by n!/(x!(n-x)!).
%
% Algorithm:
%--------------------------------------------------------------------------
% For vary small n, nCx can be computed naively as the ratio of
% factorials, using gamma(n+1) to return n!. For moderately sized n, n!
% (& x! &/or (n-x)!) become very large, and naive computation isn't
% possible. Direct computation is always possible upon noting that the
% expression cancels down to the product of x fractions:
%
%             n!       n   n-1         n-(x-1)
%         ---------  = - x --- x ... x -------
%         n! (n-x)!    x   x-1            1
%
% Unfortunately this cunning computation (given at the end of this
% function) is difficult to vectorise. Therefore we compute via the log
% of nCx as e^(ln(n!)-ln(x!)-ln((n-x)!), using the special function
% gammaln:
%
%       nCx = exp( gammaln(n+1) - gammaln(x+1) - gammaln(n-x+1) )
%
% The result is rounded to cope with roundoff error for smaller values
% of n & x. See Press et al., Sec6.1 for further details.
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
% Copyright (C) 1999-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_Ipdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<3, p=0.5; end
if nargin<2, error('Insufficient arguments'), end
ad = [ndims(x);ndims(n);ndims(p)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(n),ones(1,rd-ad(2))];...
      [size(p),ones(1,rd-ad(3))]];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size');
end

%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
f = zeros(rs);

%-Only defined for whole n, and for p in [0,1]. Return NaN if undefined.
md = ( ones(size(x))  &  n==floor(n)  &  n>=0  &  p>=0  &  p<=1 );
if any(~md(:))
    f(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Non-zero only where defined and x is whole with 0<=x<=n
Q  = find( md  &  x==floor(x)  &  n>=x  &  x>=0 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qn=Q; else Qn=1; end
if xa(3), Qp=Q; else Qp=1; end

%-Compute
f(Q) = round(exp(gammaln(n(Qn)+1) -gammaln(x(Qx)+1) - gammaln(n(Qn)-x(Qx)+1)));
f(Q) = f(Q).* p(Qp).^x(Qx) .* (1-p(Qp)).^(n(Qn)-x(Qx));

%-Return
%--------------------------------------------------------------------------
return



%==========================================================================
%-Direct computation method: (For interest)
%==========================================================================
% The following cunning direct computation is faster than using log
% gammas, but is rather difficult to vectorise.
%q=1-p;
%if r<n/2
%   %-better for small r (less terms) / small p (smaller numbers)
%   %----------------------------------------------------------------------
%   f=prod([[n:-1:n-r+1]*p,1]./[r:-1:1,1])*q^(n-r);
%else
%   %-better for large r (less terms) / small q (smaller numbers)
%   %----------------------------------------------------------------------
%   f=prod([[n:-1:r+1]*q,1]./[n-r:-1:1,1])*p^r;
%end
