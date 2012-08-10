function F = spm_Ncdf_jdw(x,u,v)
% Cumulative Distribution Function (CDF) for univariate Normal distributions: J.D.  Williams aproximation
% FORMAT F = spm_Ncdf_jdw(x,u,v)
%
% x - ordinates
% u - mean              [Defaults to 0]
% v - variance  (v>0)   [Defaults to 1]
% F - pdf of N(u,v) at x (Lower tail probability)
%__________________________________________________________________________
%
% spm_Ncdf implements the Cumulative Distribution Function (CDF) for
% the Normal (Gaussian) family of distributions.
%
% References:
%--------------------------------------------------------------------------
% An Approximation to the Probability Integral
% J. D. Williams 
% The Annals of Mathematical Statistics, Vol. 17, No. 3. (Sep., 1946), pp.
% 363-365. 
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_Ncdf_jdw.m 4836 2012-08-10 15:55:21Z karl $


%-Format arguments
%--------------------------------------------------------------------------
if nargin < 3, v = 1; end
if nargin < 2, u = 0; end

%-Approximate integral
%--------------------------------------------------------------------------
x    = (x - u)./sqrt(abs(v));
F    = sqrt(1 - exp(-(2/pi)*x.^2))/2;
i    = x < 0;
F(i) = -F(i);
F    = F + 1/2;