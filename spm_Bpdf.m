function f = spm_Bpdf(x,v,w)
% Probability Density Function (PDF) of Beta distribution
% FORMAT f = spm_Bpdf(x,v,w)
%
% x   - Beta variates (Beta has range [0,1])
% v   - Shape parameter (v>0)
% w   - Shape parameter (w>0)
% F   - PDF of Beta distribution with shape parameters [v,w] at points x
%__________________________________________________________________________
%
% spm_Bpdf implements the Probability Density Function for Beta distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The PDF of the Beta distribution shape parameters v & w, defined
% for positive integer degrees of freedom v>0 & w>0, and for x in
% [0,1] is given by: (See Evans et al., Ch5)
%
%            x^(v-1) * (1-x)^(w-1)
%    f(x) = -----------------------
%                 beta(v,w)
%
% Variate relationships:
%--------------------------------------------------------------------------
% Many: See Evans et al., Ch5
%
% Algorithm:
%--------------------------------------------------------------------------
% Direct computation using logs and MATLAB's implementation of the log
% beta function (betaln).
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
% $Id: spm_Bpdf.m 4182 2011-02-01 12:29:09Z guillaume $


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
f = zeros(rs);

%-Only defined for x in [0,1] & strictly positive v & w.
% Return NaN if undefined.
md = ( x>=0  &  x<=1  &  v>0  &  w>0 );
if any(~md(:))
    f(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Special cases: x=0 & x=1
f(md & x==0 & v<1) = Inf;
f(md & x==1 & w<1) = Inf;

%-Non-zero where defined & x in (0,1)
Q  = find( md  &  x>0  &  x<1 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end
if xa(3), Qw=Q; else Qw=1; end

%-Compute
f(Q) = exp((v(Qv)-1).*log(x(Qx)) + (w(Qw)-1).*log(1-x(Qx)) -betaln(v(Qv),w(Qw)));
