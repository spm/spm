function F = spm_Icdf(x,n,p)
% Cumulative Distribution Function (CDF) of Binomial Bin(n,p) distribution
% FORMAT F = spm_Icdf(x,n,p)
%
% x - ordinate
% n - Binomial n
% p - Binomial p [Defaults to 0.5]
% F - CDF
%__________________________________________________________________________
%
% spm_Icdf returns the Cumulative Distribution Function for the
% Binomial family of distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The Bin(n,p) distribution is the distribution of the number of
% successes from n identical independent Bernoulli trials each with
% success probability p. If random variable X is the number of
% successes from such a set of Bernoulli trials, then the CDF F(x) is
% Pr{X<=x}, the probability of x or less sucesses.
%
% The Binomial CDF is defined for whole n (i.e. non-negative integer n)
% and p in [0,1], given by: (See Evans et al., Ch6)
%
%           { 0                                         for x<0
%           |    _ floor(x)
%    F(x) = |    >      nCi * p^i * (1-p)^(n-i)         for 0<=x<=n
%           |    - i=0
%           { 1                                         for x>n
%
% where nCx is the Binomial coefficient "n-choose-x", given by n!/(x!(n-x)!)
%
% Normal approximation:
%--------------------------------------------------------------------------
% For (npq>5 & 0.1<=p<=0.9) | min(np,nq)>10 | npq>25 the Normal
% approximation to the Binomial may be used:
%       X~Bin(n,p),  X~:~N(np,npq)              ( ~:~ -> approx. distributed as)
% where q=1-p. With continuity correction this gives:
%       F(x) \approx \Phi((x+0.5-n*p)/sqrt(n*p*q))
% for Phi the standard normal CDF, related to the error function by
%       \Phi(x) = 0.5+0.5*erf(x/sqrt(2))
%
% Algorithm:
%--------------------------------------------------------------------------
% F(x), the CDF of the Binomial distribution, for X~Bin(n,p), is related
% to the incomplete beta function, by:
%
%    F(x) = 1 - betainc(p,x+1,n-x) (0<=x<n)
%
% See Abramowitz & Stegun, 26.5.24; and Press et al., Sec6.1 for
% further details.
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
% $Id: spm_Icdf.m 4182 2011-02-01 12:29:09Z guillaume $



%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<3, p=0.5; end
if nargin<2, error('Insufficient arguments'), end
ad = [ndims(x);ndims(n);ndims(p)];
rd = max(ad);
as = [  [size(x),ones(1,rd-ad(1))];...
    [size(n),ones(1,rd-ad(2))];...
    [size(p),ones(1,rd-ad(3))]     ];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size');
end

%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
F = zeros(rs);

%-Only defined for whole n, and for p in [0,1]. Return NaN if undefined.
md = ( ones(size(x))  &  n==floor(n)  &  n>=0  &  p>=0  &  p<=1 );
if any(~md(:))
    F(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-F is 1 where x>=n, or (p=0 & x>=0) (where betainc involves log of zero)
m1 = ( x>=n  |  (p==0 & x>=0) );
F(md&m1) = 1;

%-Non-zero only where defined & x>=0 & ~(p==1 & x<n)
% (Leave those already set to 1)
Q  = find( md  &  ~m1  &  x>=0  &  ~(p==1 & x<n) );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qn=Q; else Qn=1; end
if xa(3), Qp=Q; else Qp=1; end

%-F(x) is a step function with steps at {0,1,...,n}, so floor(x)
fxQx = floor(x(Qx));

%-Compute
F(Q) = 1 - betainc(p(Qp),fxQx+1,n(Qn)-fxQx);
