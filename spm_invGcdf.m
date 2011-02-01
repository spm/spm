function x = spm_invGcdf(F,h,l,tol)
% Inverse Cumulative Distribution Function (CDF) of Gamma distribution
% FORMAT x = spm_invGcdf(p,h,l,tol)
%
% F   - CDF (lower tail p-value)
% h   - Shape parameter (h>0)
% l   - Scale parameter (l>0)
% x   - Gamma ordinates at which CDF F(x)=F
% tol - tolerance [default 10^-6]
%__________________________________________________________________________
%
% spm_invGcdf implements the inverse Cumulative Distribution Function
% for Gamma distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The Gamma distribution has shape parameter h and scale l, and is
% defined for h>0 & l>0 and for x in [0,Inf) (See Evans et al., Ch18,
% but note that this reference uses the alternative parameterisation of
% the Gamma with scale parameter c=1/l)
%
% The Cumulative Distribution Function (CDF) F(x) is the probability
% that a realisation of a Gamma random variable X has value less than
% x. F(x)=Pr{X<x}:  (See Evans et al., Ch18)
%
% Variate relationships: (Evans et al., Ch8 & Ch18)
%--------------------------------------------------------------------------
% For natural (strictly +ve integer) shape h this is an Erlang distribution.
%
% The Standard Gamma distribution has a single parameter, the shape h.
% The scale taken as l=1.
%
% The Chi-squared distribution with v degrees of freedom is equivalent
% to the Gamma distribution with scale parameter 1/2 and shape parameter v/2.
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
% Copyright (C) 1992-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_invGcdf.m 4182 2011-02-01 12:29:09Z guillaume $



%-Parameters
%--------------------------------------------------------------------------
Dtol   = 10^-8;
maxIt = 10000;

%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<4, tol = Dtol; end
if nargin<3, error('Insufficient arguments'), end

ad = [ndims(F);ndims(h);ndims(l)];
rd = max(ad);
as = [[size(F),ones(1,rd-ad(1))];...
      [size(h),ones(1,rd-ad(2))];...
      [size(l),ones(1,rd-ad(3))]];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size');
end

%-Computation - Undefined and special cases
%--------------------------------------------------------------------------
%-Initialise result to zeros
x = zeros(rs);

%-Only defined for F in [0,1] & strictly positive h & l.
% Return NaN if undefined.
md = ( F>=0  &  F<=1  &  h>0  &  l>0 );
if any(~md(:))
    x(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Special cases: x=0 when F=0, x=Inf when F=1
x(md & F==1) = Inf;

%-Compute where defined & not special case
%--------------------------------------------------------------------------
Q  = find( md  &  F>0  &  F<1 );
if isempty(Q), return, end
if xa(1), FQ=F(Q); FQ=FQ(:); else FQ=F*ones(length(Q),1); end
if xa(2), hQ=h(Q); hQ=hQ(:); else hQ=h*ones(length(Q),1); end
if xa(3), lQ=l(Q); lQ=lQ(:); else lQ=l*ones(length(Q),1); end
%-?Q=?Q(:) stuff is to avoid discrepant orientations of vector arguments!

%-Compute initial interval estimates [0,mean] & Grow to encompass F(QF);
%--------------------------------------------------------------------------
a  = zeros(length(Q),1);
b  = hQ./(2*lQ);
QQ = 1:length(Q);
while ~isempty(QQ)
    b(QQ) = 2*b(QQ);
    QQ    = find(spm_Gcdf(b,hQ,lQ) < FQ);
end

%-Interval bisection
%--------------------------------------------------------------------------
fa = spm_Gcdf(a,hQ,lQ)-FQ;
fb = spm_Gcdf(b,hQ,lQ)-FQ;
i  = 0;
xQ = a+(b-a)/2;
QQ = 1:length(Q);

while ~isempty(QQ) &&  i<maxIt
    fxQQ        = spm_Gcdf(xQ(QQ),hQ(QQ),lQ(QQ))-FQ(QQ);
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
