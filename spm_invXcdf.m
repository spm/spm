function x = spm_invXcdf(F,v)
% Inverse Cumulative Distribution Function (CDF) of Chi-squared distribution
% FORMAT x = spm_invXcdf(F,v)
%
% F - CDF (lower tail p-value)
% v - degrees of freedom (v>0, non-integer d.f. accepted)
% x - Chi-squared ordinates at which CDF F(x)=F
%__________________________________________________________________________
%
% spm_invXcdf implements the inverse Cumulative Distribution of the
% Chi-squared distribution.
%
% Definition:
%--------------------------------------------------------------------------
% The Chi-squared distribution with v degrees of freedom is defined for
% positive integer v and x in [0,Inf). The Cumulative Distribution
% Function (CDF) F(x) is the probability that a realisation of a
% Chi-squared random variable X has value less than x. F(x)=Pr{X<x}:
% (See Evans et al., Ch8)
%
% Variate relationships: (Evans et al., Ch8 & Ch18)
%--------------------------------------------------------------------------
% The Chi-squared distribution with v degrees of freedom is equivalent to
% the Gamma distribution with scale parameter 1/2 and shape parameter v/2.
%
% Algorithm:
%--------------------------------------------------------------------------
% Using routine spm_invGcdf for Gamma distribution, with appropriate
% parameters.
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
% Copyright (C) 1993-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_invXcdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Check enough arguments
%--------------------------------------------------------------------------
if nargin<2, error('Insufficient arguments'), end

%-Computation
%--------------------------------------------------------------------------
x = spm_invGcdf(F,v/2,1/2);
