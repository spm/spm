function [p] = spm_dirichlet(x,alpha);
% Dirichlet distribution
% 
% FORMAT [p] = dirichlet(x,alpha)
% 
% x     - vector of outcome/event probabilities
% alpha - vector of observed events
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% 
% Will Penny
% $Id: spm_dirichlet.m 2507 2008-11-30 14:45:22Z klaas $

alpha = alpha(:);
x     = x(:);

logp  = gammaln(sum(alpha))+sum((alpha-1).*log(x+eps))-sum(gammaln(alpha));
p     = exp(logp);

return
