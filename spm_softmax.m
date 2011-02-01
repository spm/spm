function [y] = spm_softmax(x,k)
% softmax (neural transfer)  function
% FORMAT [y] = spm_softmax(x,k)
%
% x - vector of activity
% k - temperature or inverse sensitivity parameter (default k = 1)
%
% y   = exp(k*x)/sum(exp(k*x))
 
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_softmax.m 4187 2011-02-01 20:13:57Z karl $
 
% apply
%--------------------------------------------------------------------------
if nargin == 1
    k = 1; 
end
x   = x - max(x);
y   = exp(k*x)/sum(exp(k*x));

