function f = spm_Tpdf(x,v)
% Probability Density Function (PDF) of Students t distribution
% FORMAT f = spm_Tpdf(x,v)
%
% x - t-ordinates
% v - degrees of freedom (v>0, non-integer d.f. accepted)
% f - PDF of t-distribution with v degrees of freedom (df) at points t
%__________________________________________________________________________
%
% spm_Tpdf implements the Probability Density Function of Students 
% t-distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The Student's t-distribution with v degrees of freedom is defined for
% positive integer v and x in (-Inf,Inf), and has Probability Distribution
% Function (PDF) f(x) given by: (See Evans et al., Ch37)
%
%               gamma((v+1)/2)
%    f(x) = -----------------------  *  (1 + x^2/v) ^ -((v+1)/2
%           sqrt(pi*v) * gamma(v/2)
%
% This implementation is not restricted to whole (positive integer) df
% v, rather it will compute for any df v>0.
%
% Algorithm:
%--------------------------------------------------------------------------
% Direct computation using the beta function for
%       sqrt(pi)*gamma(v/2) / gamma((v+1)/2)  =  beta(v/2,1/2)
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
% $Id: spm_Tpdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<2, error('Insufficient arguments'), end

ad = [ndims(x);ndims(v)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(v),ones(1,rd-ad(2))];];
rs = max(as);
xa = prod(as,2)>1;
if all(xa) && any(diff(as(xa,:)))
    error('non-scalar args must match in size');
end


%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
f = zeros(rs);

%-Only defined for v>0. Return NaN if undefined.
md = ( ones(size(x))  &  v>0 );
if any(~md(:))
    f(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Compute where defined
Q  = find( md );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end

%-Compute
f(Q) = ( (1+x(Qx).^2./v(Qv)).^(-(v(Qv)+1)/2) ) ./ ...
     (sqrt(v(Qv)).*beta(v(Qv)/2,1/2));
