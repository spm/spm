function F = spm_Fcdf(x,v,w)
% Cumulative Distribution Function (CDF) of F (Fisher-Snedecor) distribution
% FORMAT F = spm_Fcdf(x,df)
% FORMAT F = spm_Fcdf(x,v,w)
%
% x  - F-variate   (F has range [0,Inf) )
% df - Degrees of freedom, concatenated along last dimension
%      Eg. Scalar (or column vector) v & w. Then df=[v,w];
% v  - Shape parameter 1 /   numerator degrees of freedom (v>0)
% w  - Shape parameter 2 / denominator degrees of freedom (w>0)
% F  - CDF of F-distribution with [v,w] degrees of freedom at points x
%__________________________________________________________________________
%
% spm_Fcdf implements the Cumulative Distribution Function of the F-distribution.
%
% Definition:
%--------------------------------------------------------------------------
% The CDF F(x) of the F distribution with degrees of freedom v & w,
% defined for positive integer degrees of freedom v & w, is the
% probability that a realisation of an F random variable X has value
% less than x F(x)=Pr{X<x} for X~F(v,w). The F-distribution is defined
% for v>0 & w>0, and for x in [0,Inf) (See Evans et al., Ch16).
%
% Variate relationships: (Evans et al., Ch16 & 37)
%--------------------------------------------------------------------------
% The square of a Student's t variate with w degrees of freedom is
% distributed as an F-distribution with [1,w] degrees of freedom.
%
% For X an F-variate with v,w degrees of freedom, w/(w+v*X^2) has
% distributed related to a Beta random variable with shape parameters
% w/2 & v/2.
%
% Algorithm:
%--------------------------------------------------------------------------
% Using the relationship with the Beta distribution: The CDF of the
% F-distribution with v,w degrees of freedom is related to the
% incomplete beta function by:
%       Pr(X<x) = 1 - betainc(w/(w+v*x^2),w/2,v/2)
% See Abramowitz & Stegun, 26.6.2; Press et al., Sec6.4 for
% definitions of the incomplete beta function. The relationship is
% easily verified by substituting for w/(w+v*x^2) in the integral of the
% incomplete beta function.
%
% MATLAB's implementation of the incomplete beta function is used.
%
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
% Copyright (C) 1992-2013 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_Fcdf.m 5602 2013-08-12 13:35:52Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<2, error('Insufficient arguments'), end

%-Unpack degrees of freedom v & w from single df parameter (v)
if nargin<3
    vs = size(v);
    if prod(vs)==2
        %-DF is a 2-vector
        w = v(2); v = v(1);
    elseif vs(end)==2
        %-DF has last dimension 2 - unpack v & w
        nv = prod(vs);
        w  = reshape(v(nv/2+1:nv),vs(1:end-1));
        v  = reshape(v(1:nv/2)   ,vs(1:end-1));
    else
        error('Can''t unpack both df components from single argument')
    end
end

%-Check argument sizes
ad = [ndims(x);ndims(v);ndims(w)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(v),ones(1,rd-ad(2))];...
      [size(w),ones(1,rd-ad(3))]];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size'), end

%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
F = zeros(rs);

%-Only defined for strictly positive v & w. Return NaN if undefined.
md = ( ones(size(x))  &  v>0  &  w>0 );
if any(~md(:))
    F(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Non-zero where defined and x>0
Q  = find( md  &  x>0 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end
if xa(3), Qw=Q; else Qw=1; end

%-Compute
F(Q) = 1 - betainc(w(Qw)./(w(Qw) + v(Qv).*x(Qx)),w(Qw)/2,v(Qv)/2);
