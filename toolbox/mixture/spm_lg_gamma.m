function [lng] = spm_lg_gamma (p,b)
% Log of generalised gamma function
% FORMAT [lng] = spm_lg_gamma(p,b)
%
% p       - dimension parameter
% b       - degrees of freedom type parameter
%__________________________________________________________________________
%
% References:
% * Bayesian Inference in Statistical Analysis, Box & Tiao, 1992, p. 427.
% * Muirhead p. 62.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_lg_gamma.m 2696 2009-02-05 20:29:48Z guillaume $

if b <= (p-1)/2
   warning('Parameter out of range');
   return
end

lng = (p*(p-1)/2) * gammaln(0.5);
for alpha = 1:p
  lng = lng + gammaln(b+0.5*(alpha-1));
end
