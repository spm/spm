function x = spm_invFcdf(F,v,w)
% Inverse Cumulative Distribution (CDF) of F (Fisher-Snedecor) distribution
% FORMAT x = spm_invFcdf(F,df)
% FORMAT x = spm_invFcdf(F,v,w)
%
% F  - CDF (lower tail p-value)
% df - Degrees of freedom, concatenated along last dimension
%      Eg. Scalar (or column vector) v & w. Then df=[v,w];
% v  - Shape parameter 1 /   numerator degrees of freedom (v>0)
% w  - Shape parameter 2 / denominator degrees of freedom (w>0)
% x  - F-variate   (F has range [0,Inf) )
%__________________________________________________________________________
%
% spm_Fcdf implements the inverse Cumulative Distribution Function
% for the F-distribution.
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
% distribution related to a Beta random variable with shape parameters
% w/2 & v/2, as described below.
%
% Algorithm:
%--------------------------------------------------------------------------
% Using the routine spm_invBcdf for the Beta distribution, with
% appropriate parameters:  The CDF of the F-distribution with v,w
% degrees of freedom is related to the incomplete beta function by:
%       Pr(X<x) = 1 - betainc(w/(w+v*x^2),w/2,v/2)
% See Abramowitz & Stegun, 26.6.2; Press et al., Sec6.4 for
% definitions of the incomplete beta function. The relationship is
% easily verified by substituting for w/(w+v*x^2) in the integral of the
% incomplete beta function.
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
% Copyright (C) 1993-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_invFcdf.m 6378 2015-03-15 14:46:41Z karl $


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
ad = [ndims(F);ndims(v);ndims(w)];
rd = max(ad);
as = [[size(F),ones(1,rd-ad(1))];...
      [size(v),ones(1,rd-ad(2))];...
      [size(w),ones(1,rd-ad(3))]];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size');
end

%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
x = zeros(rs);

%-Only defined for F in [0,1] & strictly positive v & w.
% Return NaN if undefined.
md = ( F>=0  &  F<=1  &  v>0  &  w>0 );
if any(~md(:))
    x(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Special cases: x=0 when F=0, x=Inf when F=1
x(md & F==1) = Inf;

%-Compute where defined & not special case
Q  = find( md  &  F>0  &  F<1 );
if isempty(Q), return, end
if xa(1), QF=Q; else QF=1; end
if xa(2), Qv=Q; else Qv=1; end
if xa(3), Qw=Q; else Qw=1; end

%-Compute
bQ   = spm_invBcdf(1-F(QF),w(Qw)/2,v(Qv)/2);
x(Q) = (w(Qw)./bQ -w(Qw))./v(Qv);
