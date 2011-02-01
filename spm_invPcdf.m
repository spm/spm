function x = spm_invPcdf(F,l)
% Inverse Cumulative Distribution Function (CDF) of Poisson distribution
% FORMAT x = spm_invPcdf(F,l)
%
% F - CDF (lower tail p-value)
% x - ordinates
% l - Poisson mean parameter (lambda l>0) [Defaults to 1]
%__________________________________________________________________________
%
% spm_invPcdf returns the inverse Cumulative Distribution Function for
% the Poisson family of distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The Poisson Po(l) distribution is the distribution of the number of
% events in unit time for a stationary Poisson process with mean
% parameter lambda=1, or equivalently rate 1/l. If random variable X is
% the number of such events, then X~Po(l), and the CDF F(x) is
% Pr({X<=x}.
%
% The Poisson distribution is discrete, defined for l in [0,Inf) and x
% in {0,1,...}, so F(x) is a discrete function. This inverse CDF
% function returns the smallest Whole x such that the F(x) equals or
% exceeds the given CDF probability F. I.e. F(x) is treated as a step
% function.
%
% Algorithm:
%--------------------------------------------------------------------------
% x is found by direct summation of the Poisson PDFs until F is exceeded.
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
% Copyright (C) 1996-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_invPcdf.m 4182 2011-02-01 12:29:09Z guillaume $



%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<2, l=1; end
if nargin<1, error('Insufficient arguments'), end

ad = [ndims(F);ndims(l)];
rd = max(ad);
as = [[size(F),ones(1,rd-ad(1))];...
      [size(l),ones(1,rd-ad(2))];];
rs = max(as);
xa = prod(as,2)>1;
if all(xa) && any(diff(as(xa,:)))
    error('non-scalar args must match in size');
end


%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
x = zeros(rs);

%-Only defined for F in [0,1] & l>0 & . Return NaN if undefined.
md = ( F>=0  &  F<=1  &  l>0 );
if any(~md(:))
    x(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Infinite where defined but F=1
mi = ( F==1 ); x(md&mi) = Inf;

%-Non-zero only where defined, & F>0
Q  = find( md  &  ~mi  &  F>0 );
if isempty(Q), return, end
if xa(1), QF=Q; else QF=1; end
if xa(2), Ql=Q; else Ql=1; end

%-Compute by directly summing Poisson PDF's for successive x & comparing with F
tx  = 0;
Ftx = spm_Ppdf(tx,l(Ql));
while any(F(QF)>Ftx)
    tx      = tx+1;
    i       = find(Ftx<F(QF));
    x(Q(i)) = x(Q(i)) + 1;
    Ftx     = Ftx + spm_Ppdf(tx,l(Ql));
end
