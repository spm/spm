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
% $Id: spm_softmax.m 6655 2015-12-23 20:21:27Z karl $
 
% apply
%--------------------------------------------------------------------------
if nargin > 1, x = k*x; end

n    = size(x);
if n(end) > 1
    ind    = cell(size(n));
    ind(:) = {':'};
    for i  = 1:n(end)
        sub       = ind;
        sub(end)  = {i};
        y(sub{:}) = spm_softmax(x(sub{:}));
    end
else
    x   = x - max(x);
    ex  = exp(x);
    y   = ex/sum(ex);
end

