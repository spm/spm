function F = spm_Ncdf(x,u,v)
% Cumulative Distribution Function (CDF) for univariate Normal distributions
% FORMAT F = spm_Ncdf(x,u,v)
%
% x - ordinates
% u - mean              [Defaults to 0]
% v - variance  (v>0)   [Defaults to 1]
% F - pdf of N(u,v) at x (Lower tail probability)
%__________________________________________________________________________
%
% spm_Ncdf implements the Cumulative Distribution Function (CDF) for
% the Normal (Gaussian) family of distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The CDF F(x) of a Normal distribution with mean u and variance v is
% the probability that a random realisation X from this distribution
% will be less than x. F(x)=Pr(X<=x) for X~N(u,v). See Evans et al.,
% Ch29 for further definitions and variate relationships.
%
% If X~N(u,v), then Z=(Z-u)/sqrt(v) has a standard normal distribution,
% Z~N(0,1). The CDF of the standard normal distribution is known as \Phi(z).
%
% (KWorsley) For extreme variates with abs(z)>6 where z=(x-u)/sqrt(v), the
% approximation \Phi(z) \approx exp(-z^2/2)/(z*sqrt(2*pi)) may be useful.
%
% Algorithm:
%--------------------------------------------------------------------------
% The CDF for a standard N(0,1) Normal distribution, \Phi(z), is
% related to the error function by: (Abramowitz & Stegun, 26.2.29)
%
%       \Phi(z) = 0.5 + erf(z/sqrt(2))/2
%
% MATLAB's implementation of the error function is used for computation.
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
% Copyright (C) 1995-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_Ncdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<3, v=1; end
if nargin<2, u=0; end
if nargin<1, F=[]; return, end
ad = [ndims(x);ndims(u);ndims(v)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(u),ones(1,rd-ad(2))];...
      [size(v),ones(1,rd-ad(3))]];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size');
end

%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
F = zeros(rs);

%-Only defined for strictly positive variance v. Return NaN if undefined.
md = ( ones(size(x))  &  ones(size(u))  &  v>0 );
if any(~md(:))
    F(~md) = NaN;
    warning('SPM:negativeVariance','Returning NaN for out of range arguments.');
end

%-Non-zero where defined
Q  = find( md );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qu=Q; else Qu=1; end
if xa(3), Qv=Q; else Qv=1; end

%-Compute
F(Q) = 0.5 + 0.5*erf((x(Qx)-u(Qu))./sqrt(2*v(Qv)));
