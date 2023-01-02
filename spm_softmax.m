function [y] = spm_softmax(x,k)
% Softmax (e.g., neural transfer) function over columns
% FORMAT [y] = spm_softmax(x,k)
%
% x - numeric array array
% k - precision, sensitivity or inverse temperature (default k = 1)
%
% y  = exp(k*x)/sum(exp(k*x))
%
% NB: If supplied with a matrix this routine will return the softmax
% function over columns - so that spm_softmax([x1,x2,..]) = [1,1,...]
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2010-2022 Wellcome Centre for Human Neuroimaging


% apply
%--------------------------------------------------------------------------
if nargin > 1,    x = k*x; end
if size(x,1) < 2; y = ones(size(x)); return, end

% exponentiate and normalise
%--------------------------------------------------------------------------
x  = exp(minus(x,max(x)));
y  = rdivide(x,sum(x));
