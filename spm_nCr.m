function c = spm_nCr(n,r)
% Combinatorics: n choose r, nCr
% FORMAT c = spm_nCr(n,r)
%
% n - Number of objects
% r - Number of objects to choose
% c - n choose r
%_______________________________________________________________________
%
% spm_nCr returns the number of ways of choosing r objects from a pool
% of n objects, without replacement, order unimportant. Equivalently:
% the number of ways of putting r objects into n indistinguishable urns
% with exclusion. These are the Binomial coefficients of Pascal's
% triangle:
%            ( n )      n!
%            |   | = ---------
%            ( r )   r! (n-r)!
%
% n & r must be whole numbers, with n>=r. Non-integer or out-of-range
% arguments return zero as nCr. Non-scalar n & r must have the same
% dimensons.
%
% Algorithm:
%-----------------------------------------------------------------------
% For vary small n, nCr can be computed naively as the ratio of
% factorials, using gamma(x+1) to return x!. For moderately sized n, n!
% (& r! &/or (n-r)!) become very large, and naive computation isn't
% possible. Direct computation is still possible upon noting that the
% expression cancels down to the product of r fractions:
%
%             n!       n   n-1         n-(r-1)
%         ---------  = - x --- x ... x -------
%         n! (n-r)!    r   r-1            1
%
% Unfortunately this cunning computation (given at the end of this
% function) is difficult to vectorise. Therefore we compute the log of
% nCr as (ln(n!)-ln(r!)-ln((n-r)!), using the log-gamma special
% function gammaln:
%
%   nCr = exp( gammaln(n+1) - gammaln(r+1) - gammaln(n-r+1) )
%
% The result is rounded to cope with roundoff error for smaller values
% of n & r. See Press et al., Sec6.1 for further details.
%
% References
%-----------------------------------------------------------------------
% Evans M, Hastings N, Peacock B (1993)
%   "Statistical Distributions"
%    2nd Ed. Wiley, New York
%
% Abramowitz M, Stegun IA, (1964)
%   "Handbook of Mathematical Functions"
%    US Government Printing Office
%
% Press WH, Teukolsky SA, Vetterling AT, Flannery BP (1992)
%   "Numerical Recipes in C"
%    Cambridge
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_nCr.m 1143 2008-02-07 19:33:33Z spm $



%-Format arguments, note & check sizes
%-----------------------------------------------------------------------
ad = [ndims(n);ndims(r)];
dc = max(ad);
as = [  [size(n),ones(1,dc-ad(1))];...
    [size(r),ones(1,dc-ad(2))]     ];
sc = max(as);
xa = prod(as,2)>1;
if (sum(xa)==2 & any(diff(as))), error('non-scalar args must match in size'), end

%-Computation
%-----------------------------------------------------------------------
%-Initialise result to zeros
c = zeros(sc);

%-Non zero only where n & r are whole and 0<=r<=n
Q = find( n==floor(n)  &  r==floor(r)  &  n>=r  &  r>=0 );
if isempty(Q), return, end
if xa(1), Qn=Q; else Qn=1; end
if xa(2), Qr=Q; else Qr=1; end

%-Compute
c(Q) = round(exp(gammaln(n(Qn)+1) -gammaln(r(Qr)+1) - gammaln(n(Qn)-r(Qr)+1)));

%-Return
%-----------------------------------------------------------------------
return


%=======================================================================
%-Direct computation method: (For interest)
%=======================================================================
% The following cunning direct computation is faster than using log
% gammas, but is rather difficult to vectorise.
%if r<n/2
%   %-better for small r (less terms), append 1's for 0!
%   c = prod([n:-1:n-r+1,1]./[r:-1:1,1]);
%else   
%   %-better for large r (less terms), append 1's for 0!
%   c = prod([n:-1:r+1,1]./[n-r:-1:1,1]);
%end
