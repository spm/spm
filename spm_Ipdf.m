function f=spm_Bpdf(r,n,p);
% Probability Distribution Function of Binomial distribution
% FORMAT f=spm_Bpdf(r,n,p);
%
% r - ordinate
% n - Binomial n
% p - Binomial p [Defaults to 0.5]
% f - PDF
%_______________________________________________________________________
%
% spm_Bpdf returns the Probability (Distribution) Function (PDF) for
% the Binomial family of distributions.
%
% Accurate computation avoiding rounding errors and large factorials is
% achieved by cunning computation.
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%-Condition & check arguments
%-----------------------------------------------------------------------
if nargin<3, p=0.5; end
if nargin<2, error('Insufficient arguments'); end
if any([size(r),size(n),size(p)]>1), error('Can''t handle vector r/n/p'), end
if any(floor([r,n])~=ceil([r,n])), error('r/n must be integers'), end

%-Out of range values
%-----------------------------------------------------------------------
if r<0 | r>n, f=0; return, end

%-Computation
%-----------------------------------------------------------------------
q=1-p;

if r<n/2
	%-better for small r (less terms) / small p (smaller numbers)
	%---------------------------------------------------------------
	f=prod([[n:-1:n-r+1]*p,1]./[r:-1:1,1])*q^(n-r);
else
	%-better for large r (less terms) / small q (smaller numbers)
	%---------------------------------------------------------------
	f=prod([[n:-1:r+1]*q,1]./[n-r:-1:1,1])*p^r;
end
