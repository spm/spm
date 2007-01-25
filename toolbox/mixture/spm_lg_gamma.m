function [lng] = spm_lg_gamma (p,b)
% Log of generalised gamma function. See eg. Box and Tiao p.427
% FORMAT [lng] = spm_lg_gamma(p,b)
%
% p          dimension parameter
% b          degrees of freedom type parameter
%___________________________________________________________________________
% Copyright (C) 2007 Wellcome Department of Imaging Neuroscience

% Will Penny 
% $Id$

if ~(b > 0.5*(p-1))
   disp('Warning in log_gen_gamma: parameter out of range');
   return
end

lng=(0.5*p*(p-1))*gammaln(0.5);
for alpha=1:p,
  % Muirhead (p.62):
  lng=lng+gammaln(b+0.5*(alpha-1));
end



