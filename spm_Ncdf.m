function p = spm_Ncdf(x,Mu,SigmaSq)
% Cumulative Distribution Function (CDF) for univariate normal distributions
% FORMAT p = spm_Ncdf(x,Mu,Var)
%
% x       - N(Mu,SigmaSq) variates
% Mu      - Mean of the univariate Normal distribution
% SigmaSq - Variance of the univariate Normal distribution
% p       - lower tail probability
%__________________________________________________________________________
%
% spm_Ncdf implements the Cumulative Distribution Function (CDF) for
% the Normal (Gaussian) family of distributions.
%
% Returns the probability p, that a variate from a Normal distribution with
% mean Mu and variance SigmaSq is less than z.
% p = Pr(Z <= z) for Z ~ N(Mu,SigmaSq).
%
% The CDF for a standard N(0,1) Normal distribution (known as \Phi), is
% related to the error function, and MatLab's implementation of the error
% function is used for computation. See "Numerical Recipies in C"
%
% For extreme variates with abs(z)>6 where z=(x-Mu)/sqrt(SigmaSq), the
% approximation \phi(z) \approx exp(-z^2/2)/(z*sqrt(2*pi)) is used
%
%---------------------------------------------------------------------------
% %W% Andrew Holmes %E%

%-Parameters - (Absolute) threshold above which the high value
%		to \phi approximation is used
%-----------------------------------------------------------------------
z1 = 6;

%-Condition arguments
%---------------------------------------------------------------------------
if nargin < 3, SigmaSq = 1;  end
if nargin < 2, Mu      = 0;  end
if nargin < 1, z       = []; end

if SigmaSq <=0 , error('SigmaSq must be strictly positive'), end

%-Computation
%=======================================================================
p  = zeros(size(x));

%-Center and scale to standard Gaussian variates
%-----------------------------------------------------------------------
z  = (x - Mu)/sqrt(SigmaSq);

%-Find really extreme variates
%-----------------------------------------------------------------------
h  = find(abs(z) >  z1);
l  = find(abs(z) <= z1);

%-Standard Ncdf
%-----------------------------------------------------------------------
if length(l)
	p(l) = 0.5 + 0.5*erf(z(l)/sqrt(2));
end

%-High value approximation to Ncdf
%-----------------------------------------------------------------------
if length(h)
	p(h) = (z(h)>0) - exp(-(z(h).^2)/2)./(z(h)*sqrt(2*pi));
end
