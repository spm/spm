function p = spm_Ncdf(z,Mu,SigmaSq)
% Cumulative Distribution Function (CDF) for univariate normal distributions
% FORMAT p = spm_Ncdf(z,Mu,Var)
%
% z       - N(Mu,SigmaSq) variates
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
%---------------------------------------------------------------------------
% %W% Andrew Holmes %E%

% - Andrew Holmes - V1 - 09/92 - ah_phi, for standard Normal variates
%                   V2 - 02/95 - Rewritten for SPM, including Mu and SigmaSq

%-Condition arguments
%---------------------------------------------------------------------------
if nargin<3, SigmaSq=1; end
if nargin<2, Mu=0; end
if nargin<1, z=[]; end

if SigmaSq<=0, error('SigmaSq must be strictly positive'), end

%-Computation
%---------------------------------------------------------------------------
if isempty(z), p=[]; return, end
p = 0.5 + 0.5*erf((z - Mu)./sqrt(2*SigmaSq));
