function x = spm_invTcdf(F,v)
% Inverse Cumulative Distribution Function (CDF) of Students t distribution
% FORMAT x = spm_invTcdf(F,v)
%
% F  - CDF (lower tail p-value)
% v - degrees of freedom (v>0, non-integer d.f. accepted)
% x - T-variate (Student's t has range (-Inf,Inf))
%__________________________________________________________________________
%
% spm_invTcdf implements the inverse Cumulative Distribution of Students
% t-distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The Student's t-distribution with v degrees of freedom is defined for
% positive integer v and real x. The Cumulative Distribution
% Function (CDF) F(x) is the probability that a realisation of a
% t-distributed random variable X has value less than x. F(x)=Pr{X<x}:
% (See Evans et al., Ch37)
%
% This implementation is not restricted to whole (positive integer) df
% v, rather it will compute for any df v>0.
%
% Variate relationships: (Evans et al., Ch37 & 7)
%--------------------------------------------------------------------------
% The Student's t distribution with 1 degree of freedom is the Standard
% Cauchy distribution, which has a simple closed form CDF.
%
% For X a t-variate with v degrees of freedom, v/(v+X^2) has
% distribution related to a Beta random variable with shape parameters
% w/2 & 1/2, as described below.
%
% Algorithm:
%--------------------------------------------------------------------------
% Using the routine spm_invBcdf for the Beta distribution, with
% appropriate parameters:  The CDF of the Student's t-distribution with
% v degrees of freedom is related to the incomplete beta function by:
%       Pr(|X|<x) = betainc(v/(v+x^2),v/2,1/2)
% so
%              {     betainc(v/(v+x^2),v/2,1/2) / 2      for x<0
%       F(x) = |   0.5                                   for x=0
%              { 1 - betainc(v/(v+x^2),v/2,1/2) / 2      for x>0
%
% See Abramowitz & Stegun, 26.5.27 & 26.7.1; Press et al., Sec6.4 for
% definitions of the incomplete beta function. The relationship is
% easily verified by substituting for v/(v+x^2) in the integral of the
% incomplete beta function.
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
% $Id: spm_invTcdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<2, error('Insufficient arguments'), end

ad = [ndims(F);ndims(v)];
rd = max(ad);
as = [  [size(F),ones(1,rd-ad(1))];...
    [size(v),ones(1,rd-ad(2))]     ];
rs = max(as);
xa = prod(as,2)>1;
if all(xa) && any(diff(as(xa,:)))
    error('non-scalar args must match in size');
end


%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
x = zeros(rs);

%-Only defined for F in [0,1] & strictly positive v.
% Return NaN if undefined.
md = ( F>=0  &  F<=1  &  v>0 );
if any(~md(:))
    x(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Special case: x is 0 when F=0.5, -Inf when F=0, +Inf when F=1
x(md & F==0) = -Inf;
x(md & F==1) = +Inf;

%-Special case: Standard Cauchy distribution when v=1
ml = ( md  &  v==1 ); if xa(1), mlF=ml; else mlF=1; end
x(ml) = tan(pi*(F(mlF)-0.5));

%-Compute where defined & not special cases
Q  = find( md  &  F>0  &  F~=0.5  &  F<1  &  v~=1 );
if isempty(Q), return, end
if xa(1), QF=Q; else QF=1; end
if xa(2), Qv=Q; else Qv=1; end

%-Compute
xQPos = F(QF)>0.5;
bQ    = spm_invBcdf(2*(xQPos -(xQPos*2-1).*F(QF)),v(Qv)/2,1/2);
x(Q)  = (xQPos*2-1) .* sqrt(v(Qv)./bQ -v(Qv));
