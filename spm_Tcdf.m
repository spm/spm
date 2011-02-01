function F = spm_Tcdf(x,v)
% Cumulative Distribution Function (CDF) of Students t distribution
% FORMAT p = spm_Tcdf(x,v)
%
% x - T-variate (Student's t has range (-Inf,Inf)
% v - degrees of freedom (v>0, non-integer d.f. accepted)
% F - CDF of Student's t-distribution with v degrees of freedom at points x
%__________________________________________________________________________
%
% spm_Tcdf implements the Cumulative Distribution of the Students t-distribution.
%
% Definition:
%--------------------------------------------------------------------------
% The CDF F(x) of the Student's t-distribution with v degrees of
% freedom is the probability that a realisation of a t random variable
% X has value less than x; F(x)=Pr{X<x} for X~G(h,c). Student's
% t-distribution is defined for real x and positive integer v (See
% Evans et al., Ch37).
%
% This implementation is not restricted to whole (positive integer) df
% v, rather it will compute for any df v>0.
%
% Variate relationships: (Evans et al., Ch37 & 7)
%--------------------------------------------------------------------------
% The Student's t distribution with 1 degree of freedom is the Standard
% Cauchy distribution, which has a simple closed form CDF.
%
% Algorithm:
%--------------------------------------------------------------------------
% The CDF of the Student's t-distribution with v degrees of freedom
% is related to the incomplete beta function by:
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
% Copyright (C) 1992-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_Tcdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<2, error('Insufficient arguments'), end

ad = [ndims(x);ndims(v)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(v),ones(1,rd-ad(2))]];
rs = max(as);
xa = prod(as,2)>1;
if all(xa) && any(diff(as(xa,:)))
    error('non-scalar args must match in size');
end


%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
F = zeros(rs);

%-Only defined for strictly positive v. Return NaN if undefined.
md = ( ones(size(x))  &  v>0 );
if any(~md(:))
    F(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Special case: f is 0.5 when x=0 (where betainc involves log of zero)
F( md  &  x==0 ) = 0.5;

%-Special case: Standard Cauchy distribution when v=1
ml = ( md  &  v==1 ); if xa(1), mlx=ml; else mlx=1; end
F(ml) = 0.5 + atan(x(mlx))/pi;

%-Compute where defined & not special cases
Q  = find( md  &  x~=0  &  v~=1 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end

%-Compute
xQxPos = x(Qx)>0;
F(Q) = xQxPos -(xQxPos*2-1).*0.5.*betainc(v(Qv)./(v(Qv)+x(Qx).^2),v(Qv)/2,1/2);
