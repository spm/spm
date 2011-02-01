function f = spm_Npdf(x,u,v)
% Probability Density Function (PDF) of univariate Normal distribution
% FORMAT f = spm_Npdf(x,u,v)
%
% x - ordinates
% u - mean              [Defaults to 0]
% v - variance  (v>0)   [Defaults to 1]
% f - pdf of N(u,v) at x
%__________________________________________________________________________
%
% spm_Npdf returns the Probability Density Function (PDF) for the
% univariate Normal (Gaussian) family of distributions.
%
% Definition:
%--------------------------------------------------------------------------
% Let random variable X have a Normal distribution with mean u and
% variance v, then Z~N(u,v). The Probability Density Function (PDF) of
% the Normal (sometimes called Gaussian) family is f(x), defined on all
% real x, given by: (See Evans et al., Ch29)
%
%                 1           ( (x-u)^2 )
%    f(r) = ------------ x exp| ------  |
%           sqrt(v*2*pi)      (   2v    )
%
% The PDF of the standard Normal distribution, with zero mean and unit
% variance, Z~N(0,1), is commonly referred to as \phi(z).
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
% Copyright (C) 1994-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_Npdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<3, v=1; end
if nargin<2, u=0; end
if nargin<1, f=[]; return, end
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
f = zeros(rs);

%-Only defined for strictly positive variance v. Return NaN if undefined.
md = ( ones(size(x))  &  ones(size(u))  &  v>0 );
if any(~md(:))
    f(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Non-zero where defined
Q  = find( md );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qu=Q; else Qu=1; end
if xa(3), Qv=Q; else Qv=1; end

%-Compute
f(Q) = exp( -(x(Qx)-u(Qu)).^2 ./ (2*v(Qv)) ) ./ sqrt(2*pi*v(Qv));
