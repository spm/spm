function z = spm_invNcdf(p,Mu,SigmaSq)
% Inverse CDF for univariate normal distributions
% FORMAT z = spm_invNcdf(p,Mu,SigmaSq)
%
% p       - lower tail probability
% Mu      - Mean of the univariate Normal distribution [defaults to 0]
% SigmaSq - Variance of the univariate Normal distribution [defaults to 1]
% z       - N(Mu,SigmaSq) variates
%_______________________________________________________________________
%
% ah_invNcdf implements the inverse of the Cumulative Distribution
% Function (CDF) for the Normal (Gaussian) family of distributions.
%
% Returns the variate z, such that Pr(Z <= z) = p for Z ~ N(Mu,SigmaSq), a
% univariate random variable distributed Normally with mean Mu and
% variance SigmaSq.
%
% The inverse CDF for a standard N(0,1) Normal distribution (known as
% \Phi^{-1}), is related to the inverse error function, and MatLab's
% implementation of the inverse error function is used for computation.
% See "Numerical Recipies in C"
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%


%-Condition arguments
%-----------------------------------------------------------------------
if nargin<3, SigmaSq=1; end
if nargin<2, Mu=0; end
if nargin<1, p=[]; end

if any( abs(p(:)-0.5) > 0.5), error('p must be in [0,1]!'), end

%-Computation
%=======================================================================
if isempty(p), z=[]; return, end
z = ( sqrt(2)*erfinv(2*p-1) .* SigmaSq ) + Mu;
