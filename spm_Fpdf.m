function f = spm_Fpdf(x,v,w)
% Probability Density Function (PDF) of F (Fisher-Snedecor) distribution
% FORMAT f = spm_Fpdf(x,df)
% FORMAT f = spm_Fpdf(x,v,w)
%
% x  - F-variate   (F has range [0,Inf) )
% df - Degrees of freedom, concatenated along last dimension
%      Eg. Scalar (or column vector) v & w. Then df=[v,w];
% v  - Shape parameter 1 /   numerator degrees of freedom (v>0)
% w  - Shape parameter 2 / denominator degrees of freedom (w>0)
% f  - PDF of F-distribution with [v,w] degrees of freedom at points x
%__________________________________________________________________________
%
% spm_Fpdf implements the Probability Density Function of the F-distribution.
%
% Definition:
%--------------------------------------------------------------------------
% The PDF of the F-distribution with degrees of freedom v & w, defined
% for positive integer degrees of freedom v>0 & w>0, and for x in
% [0,Inf) by: (See Evans et al., Ch16)
%
%             gamma((v+w)/2)  * (v/w)^(v/2) x^(v/2-1)
%    f(x) = --------------------------------------------
%           gamma(v/2)*gamma(w/2) * (1+(v/w)x)^((v+w)/2)
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
% Direct computation using the beta function for
%       gamma(v/2)*gamma(w/2) / gamma((v+w)/2)  =  beta(v/2,w/2)
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
% $Id: spm_Fpdf.m 4182 2011-02-01 12:29:09Z guillaume $


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
    error('non-scalar args must match in size');
end

%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
f = zeros(rs);

%-Only defined for strictly positive v & w. Return NaN if undefined.
md = ( ones(size(x))  &  v>0  &  w>0 );
if any(~md(:))
    f(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Non-zero where defined and x>0
Q  = find( md  &  x>0 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end
if xa(3), Qw=Q; else Qw=1; end

%-Compute
f(Q) = (v(Qv)./w(Qw)).^(v(Qv)/2) .* x(Qx).^(v(Qv)/2-1) ./ ...
    (1+(v(Qv)./w(Qw)).*x(Qx)).^((v(Qv)+w(Qw))/2) ./ ...
        beta(v(Qv)/2,w(Qw)/2);
