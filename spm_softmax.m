function [y] = spm_softmax(x,k)
% softmax (e.g., neural transfer) function over columns
% FORMAT [y] = spm_softmax(x,k)
%
% x - numeric array array
% k - precision, sensitivity or inverse temperature (default k = 1)
%
% y  = exp(k*x)/sum(exp(k*x))
%
% NB: If supplied with a matrix this rotine will return the softmax
% function over colums - so that spm_softmax([x1,x2,..]) = [1,1,...]
 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_softmax.m 7306 2018-05-07 13:42:02Z karl $
 
% apply
%--------------------------------------------------------------------------
if nargin > 1, x = k*x; end

x  = exp(bsxfun(@minus,x,max(x)));
y  = bsxfun(@rdivide,x,sum(x));
