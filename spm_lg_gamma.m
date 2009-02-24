function [lng] = spm_lg_gamma(p,b)
% Log of generalised gamma function
% FORMAT [lng] = spm_lg_gamma(p,b)
%
% p       - dimension parameter
% b       - degrees of freedom type parameter
%__________________________________________________________________________
%
% References:
% * Bayesian Inference in Statistical Analysis, Box & Tiao, 1992, p. 427.
% * Aspects of Multivariate Statistical Theory, R.J. Muirhead, p. 62.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_lg_gamma.m 2784 2009-02-24 19:11:20Z guillaume $

if b <= (p-1)/2
   warning('Parameter out of range');
   lng = NaN;
   return
end

lng = (p*(p-1)/2) * gammaln(0.5);
for alpha = 1:p
  lng = lng + gammaln(b+0.5*(alpha-1));
end
