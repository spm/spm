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
% $Id: spm_softmax.m 4579 2011-12-02 20:21:07Z karl $
 
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
    x   = x - max(x);
    y   = exp(k*x)/sum(exp(k*x));
end

