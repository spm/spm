function f = spm_Xpdf(x,v)
% Probability Density Function (PDF) of Chi-Squared distribution
% FORMAT f = spm_Xpdf(x,v)
%
% x - Chi-squared variate
% v - degrees of freedom (v>0, non-integer d.f. accepted)
% f - PDF at x of Chi-squared distribution with v degrees of freedom
%__________________________________________________________________________
%
% spm_Xpdf implements the Probability Density Function of the 
% Chi-squared distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The Chi-squared distribution with v degrees of freedom is defined for
% positive integer v and x in [0,Inf), and has Probability Distribution
% Function (PDF) f(x) given by: (See Evans et al., Ch8)
%
%           x^((v-2)/2) exp(-x/2)
%    f(x) = ---------------------
%           2^(v/2) * gamma(v/2)
%
% Variate relationships: (Evans et al., Ch8 & Ch18)
%--------------------------------------------------------------------------
% The Chi-squared distribution with v degrees of freedom is equivalent
% to the Gamma distribution with scale parameter 1/2 and shape parameter v/2.
%
% Algorithm:
%--------------------------------------------------------------------------
% Using routine spm_Gpdf for Gamma distribution, with appropriate parameters.
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
% $Id: spm_Xpdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Check enough arguments
%--------------------------------------------------------------------------
if nargin<2, error('Insufficient arguments'), end

%-Computation
%--------------------------------------------------------------------------
f = spm_Gpdf(x,v/2,1/2);
