function [y] = spm_softmax(x,k)
% softmax (neural transfer)  function
% FORMAT [y] = spm_softmax(x)
%
% x - vector of activity
% k - temperature or inverse sensitivity parameter (default k = 1)
%
% y   = exp(k*x)/sum(exp(k*x))
 
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_softmax.m 3757 2010-03-08 11:41:53Z guillaume $
 
% apply
%--------------------------------------------------------------------------
if nargin == 1
    k = 1; 
end
x   = x - max(x);
y   = exp(k*x)/sum(exp(k*x));

