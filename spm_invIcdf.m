function r = spm_invIcdf(F,n,p)
% Inverse Cumulative Distribution Function (CDF) of Binomial distribution
% FORMAT r = spm_invIcdf(F,n,p)
%
% F - CDF (lower tail p-value)
% n - Binomial n
% p - Binomial p [Defaults to 0.5]
% r - ordinate
%__________________________________________________________________________
%
% spm_invIcdf returns the inverse Cumulative Distribution Function for
% the Binomial family of distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The Bin(n,p) distribution is the distribution of the number of
% successes from n identical independent Bernoulli trials each with
% success probability p. If random variable R is the number of
% successes from such a set of Bernoulli trials, then the CDF F(r) is
% the probability of r or less sucesses.
%
% The Binomial distribution is discrete, defined for p in [0,1] and r
% in {0,1,...,n}, so F(r) is a discrete function. This inverse CDF
% function returns the smallest Whole r such that the F(r) equals or
% exceeds the given CDF probability F. I.e. F(r) is treated as a step
% function.
%
% Algorithm:
%--------------------------------------------------------------------------
% r is found by direct summation of the Binomial PDFs until F is exceeded.
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
% $Id: spm_invIcdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<3, p=0.5; end
if nargin<2, error('Insufficient arguments'), end
ad = [ndims(F);ndims(n);ndims(p)];
rd = max(ad);
as = [[size(F),ones(1,rd-ad(1))];...
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
r = zeros(rs);

%-Only defined for whole n, p&F in [0,1]. Return NaN if undefined.
md = ( F>=0  &  F<=1  &  n==floor(n)  &  n>=0  &  p>=0  &  p<=1 );
if any(~md(:))
    r(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Non-zero only where defined, & F>0
Q  = find( md  &  F>0 );
if isempty(Q), return, end
if xa(1), QF=Q; else QF=1; end
if xa(2), Qn=Q; else Qn=1; end
if xa(3), Qp=Q; else Qp=1; end

%-Compute by directly summing Bin PDF's for successive r & comparing with F
tr  = 0;
Ftr = spm_Ipdf(tr,n(Qn),p(Qp));
while any(F(QF)>Ftr) && any(n(Qn)>tr)
    tr      = tr+1;
    i       = find(Ftr<F(QF));
    r(Q(i)) = r(Q(i)) + 1;
    Ftr     = Ftr + spm_Ipdf(tr,n(Qn),p(Qp));
end
