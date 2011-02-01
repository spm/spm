function x = spm_invNcdf(F,u,v)
% Inverse Cumulative Distribution Function (CDF) for univariate Normal
% FORMAT x = spm_invNcdf(F,u,v)
%
% F - CDF (lower tail p-value)
% u - mean              [Defaults to 0]
% v - variance  (v>0)   [Defaults to 1]
% x - ordinates of N(u,v) at which CDF F(x)=F
%__________________________________________________________________________
%
% spm_invNcdf implements the inverse of the Cumulative Distribution
% Function (CDF) for the Normal (Gaussian) family of distributions.
%
% Returns the variate x, such that Pr{X<x} = F for X~N(u,v), a
% univariate random variable distributed Normally with mean u and
% variance v.
%
% Definition:
%--------------------------------------------------------------------------
% The CDF F(x) of a Normal distribution with mean u and variance v is
% the probability that a random realisation X from this distribution
% will be less than x. F(x)=Pr(X<=x) for X~N(u,v). The inverse CDF
% returns the normal ordinate x for which the CDF is F. See Evans et
% al., Ch29 for further definitions and variate relationships.
%
% If X~N(u,v), then Z=(Z-u)/sqrt(v) has a standard normal distribution,
% Z~N(0,1). The CDF of the standard normal distribution is known as \Phi(z),
% its inverse as \Phi^{-1}(F).
%
% Algorithm:
%--------------------------------------------------------------------------
% The CDF for a standard N(0,1) Normal distribution, \Phi(z), is
% related to the error function by: (Abramowitz & Stegun, 26.2.29)
%
%       \Phi(z)      = 0.5 + erf(z/sqrt(2))/2
% so
%
%       \Phi^{-1}(p) = sqrt(2) * erfinv(2*p-1)
%
% where erfinv(.) is the inverse error function.
%
% MATLAB's implementation of the inverse error function is used for
% computation of z=\Phi^{-1}(F), the corresponding standard normal
% variate, which converted to a variate x from a N(u,v) distribution by:
%       x = u+z*sqrt(v)
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
% Copyright (C) 1992-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_invNcdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<3, v=1; end
if nargin<2, u=0; end
if nargin<1, x=[]; return, end
ad = [ndims(F);ndims(u);ndims(v)];
rd = max(ad);
as = [  [size(F),ones(1,rd-ad(1))];...
    [size(u),ones(1,rd-ad(2))];...
    [size(v),ones(1,rd-ad(3))]     ];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size');
end

%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
x = zeros(rs);

%-Only defined for F in [0,1], & strictly positive variance v.
% Return NaN if undefined.
md = ( F>=0  &  F<=1  &  ones(size(u))  &  v>0 );
if any(~md(:))
    x(~md) = NaN;
    warning('SPM:outOfRangeNormal','Returning NaN for out of range arguments');
end

%-Compute where defined
Q  = find( md );
if isempty(Q), return, end
if xa(1), QF=Q; else QF=1; end
if xa(2), Qu=Q; else Qu=1; end
if xa(3), Qv=Q; else Qv=1; end

%-Compute
x(Q) = ( sqrt(2)*erfinv(2*F(QF)-1) .* sqrt(v(Qv)) ) + u(Qu);
