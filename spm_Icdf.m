function F=spm_Bcdf(r,n,p,verbose);
% Cumulative Distribution Function (CDF) of Binomial Bin(n,p) distribution
% FORMAT F=spm_Bcdf(r,n,p);
%
% r - ordinate
% n - Binomial n
% p - Binomial p [Defaults to 0.5]
% F - CDF
%_______________________________________________________________________
%
% spm_Bcdf returns the Cumulative Distribution Function for the
% Binomial family of distributions.
%
% Accurate computation avoiding rounding errors and large factorials is
% achieved by cunning computation. For large n (&r) the normal
% approximation to the Binomial distribution is used.
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%-Condition & check arguments
%-----------------------------------------------------------------------
if nargin<4 verbose=0; end
if nargin<3, p=0.5; end
if nargin<2, error('Insufficient arguments'); end
if any([size(r),size(n),size(p)]>1), error('Can''t handle vector r/n/p'), end
if any(floor([r,n])~=ceil([r,n])), error('r/n must be integers'), end

%-Out of range values
%-----------------------------------------------------------------------
if r<0 F=0; return, end
if r>n F=1; return, end

%-Computation
%-----------------------------------------------------------------------
q=1-p;

if (n>9*q/p)&(n>100)&(r>50)
	%-Normal approximation to Binomial (have large n, & r away from 0)
	% By CLT R~Bin(n,p), R~:~N(np,npq)
	%---------------------------------------------------------------
	%-Normalise r as x, a standard normal variate
	% & add .5 continuity correction
	x=(r+0.5-n*p)/sqrt(n*p*q);
	% Use phi, standard normal CDF
	F = 0.5+0.5*erf(x/sqrt(2));
else
	%-Sum Binomial PDF
	%---------------------------------------------------------------
	F=0;
	for i=0:r, F = F + spm_Bpdf(i,n,p); end
end
