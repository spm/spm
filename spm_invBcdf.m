function x = spm_invBcdf(F,v,w,tol)
% Inverse Cumulative Distribution Function (CDF) of Beta distribution
% FORMAT x = spm_invBcdf(F,v,w,tol)
%
% F   - CDF (lower tail p-value)
% v   - Shape parameter (v>0)
% w   - Shape parameter (w>0)
% x   - Beta ordinates at which CDF F(x)=F
% tol - tolerance [default 10^-6]
%__________________________________________________________________________
%
% spm_invBcdf implements the inverse Cumulative Distribution Function
% for Beta distributions.
%
% Returns the Beta-variate x such that Pr{X<x} = F for X a Beta
% random variable with shape parameters v and w.
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
% Interval bisection to find the inverse CDF to within tol.
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
% $Id: spm_invBcdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Parameters
%--------------------------------------------------------------------------
Dtol  = 10^-8;
maxIt = 10000;

%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<4, tol = Dtol; end
if nargin<3, error('Insufficient arguments'), end

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

%-Computation - Undefined and special cases
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

%-Special cases: x=0 when F=0, x=1 when F=1
x(md & F==1) = 1;

%-Compute where defined & not special case
%--------------------------------------------------------------------------
Q  = find( md  &  F>0  &  F<1 );
if isempty(Q), return, end
if xa(1), FQ=F(Q); FQ=FQ(:); else FQ=F*ones(length(Q),1); end
if xa(2), vQ=v(Q); vQ=vQ(:); else vQ=v*ones(length(Q),1); end
if xa(3), wQ=w(Q); wQ=wQ(:); else wQ=w*ones(length(Q),1); end

%-Interval bisection
%--------------------------------------------------------------------------
a  = zeros(length(Q),1); fa = a-FQ;
b  =  ones(length(Q),1); fb = b-FQ;
i  = 0;
xQ = a+1/2;
QQ = 1:length(Q);

while ~isempty(QQ) &&  i<maxIt
    fxQQ        = betainc(xQ(QQ),vQ(QQ),wQ(QQ))-FQ(QQ);
    mQQ         = fa(QQ).*fxQQ > 0;
    a(QQ(mQQ))  = xQ(QQ(mQQ));   fa(QQ(mQQ))  = fxQQ(mQQ);
    b(QQ(~mQQ)) = xQ(QQ(~mQQ));  fb(QQ(~mQQ)) = fxQQ(~mQQ);
    xQ(QQ)      = a(QQ) + (b(QQ)-a(QQ))/2;
    QQ          = QQ( (b(QQ)-a(QQ))>tol );
    i           = i+1;
end

if i==maxIt
    warning('convergence criteria not reached - maxIt reached');
end

x(Q) = xQ;
