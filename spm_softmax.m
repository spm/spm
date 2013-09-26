function [y] = spm_softmax(x,k)
% softmax (neural transfer) function of COLUMN vectors
% FORMAT [y] = spm_softmax(x,k)
%
% x - vector of activity
% k - temperature or inverse sensitivity parameter (default k = 1)
%
% y   = exp(k*x)/sum(exp(k*x))
%
% NB: If supplied with a matric this rotine will return the softmax
% function over colums - so that spm_softmax([x1,x2,..]) = [1,1,...]
 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_softmax.m 5657 2013-09-26 16:53:40Z karl $
 
% apply
%--------------------------------------------------------------------------
if nargin == 1
    k = 1;
end
n    = size(x,2);
if n > 1
    for i = 1:n
        y(:,i) = spm_softmax(x(:,i),k);
    end
else
    x   = k*(x - max(x));
    y   = exp(x)/sum(exp(x));
end

