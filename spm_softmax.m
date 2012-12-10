function [y] = spm_softmax(x,k)
% softmax (neural transfer) function of column vectors
% FORMAT [y] = spm_softmax(x,k)
%
% x - vector of activity
% k - temperature or inverse sensitivity parameter (default k = 1)
%
% y   = exp(k*x)/sum(exp(k*x))
 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_softmax.m 5105 2012-12-10 12:31:22Z karl $
 
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

