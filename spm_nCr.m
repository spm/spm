function c=spm_nCr(n,r);
% Combinatorics: n choose r
% FORMAT c=spm_nCr(n,r);
%_______________________________________________________________________
%
% spm_nCr returns the number of ways of choosing r objects from a pool
% of n objects, without replacement, order unimportant. Equivalently
% the number of ways of putting r objects into n indistinguishable urns
% with exclusion. These are the Binomial coefficients of Pascal's
% triangle.
%
% Accurate computation of Binomial coefficients nCr avoiding rounding
% errors and large factorials by cunning computation.
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%-Condition & check arguments
%-----------------------------------------------------------------------
if any([size(n),size(r)]>1), error('Can''t handle vector n/r'), end
if any(floor([n,r])~=ceil([n,r])), error('n/r must be integers'), end

%-Out of range values
%-----------------------------------------------------------------------
if ( r<0 | r>n ), error('r out of n range'); return, end

%-Computation
%-----------------------------------------------------------------------
if r<n/2
	%-better for small r (less terms)
	% append 1's for 0!
	c = prod([n:-1:n-r+1,1]./[r:-1:1,1]);
else	
	%-better for large r (less terms)
	c = prod([n:-1:r+1,1]./[n-r:-1:1,1]);
end
