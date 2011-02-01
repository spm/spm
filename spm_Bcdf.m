function F = spm_Bcdf(x,v,w)
% Inverse Cumulative Distribution Function (CDF) of Beta distribution
% FORMAT F = spm_Bcdf(x,v,w)
%
% x   - Beta variates (Beta has range [0,1])
% v   - Shape parameter (v>0)
% w   - Shape parameter (w>0)
% F   - CDF of Beta distribution with shape parameters [v,w] at points x
%__________________________________________________________________________
%
% spm_Bcdf implements the Cumulative Distribution Function for Beta
% distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The Beta distribution has two shape parameters, v and w, and is
% defined for v>0 & w>0 and for x in [0,1] (See Evans et al., Ch5).
% The Cumulative Distribution Function (CDF) F(x) is the probability
% that a realisation of a Beta random variable X has value less than
% x. F(x)=Pr{X<x}: This function is usually known as the incomplete Beta
% function. See Abramowitz & Stegun, 26.5; Press et al., Sec6.4 for
% definitions of the incomplete beta function.
%
% Variate relationships:
%--------------------------------------------------------------------------
% Many: See Evans et al., Ch5
%
% Algorithm:
%--------------------------------------------------------------------------
% Using MATLAB's implementation of the incomplete beta finction (betainc).
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
%__________________________________________________________________________
% Copyright (C) 1999-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_Bcdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<3, error('Insufficient arguments'), end

ad = [ndims(x);ndims(v);ndims(w)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
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
F = zeros(rs);

%-Only defined for x in [0,1] & strictly positive v & w.
% Return NaN if undefined.
md = ( x>=0  &  x<=1  &  v>0  &  w>0 );
if any(~md(:))
    F(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Special cases: F=1 when x=1
F(md & x==1) = 1;

%-Non-zero where defined & x>0, avoid special cases
Q  = find( md  &  x>0  &  x<1 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end
if xa(3), Qw=Q; else Qw=1; end

%-Compute
F(Q) = betainc(x(Qx),v(Qv),w(Qw));
