function [X,s] = spm_ndgrid(x)
% returns a matrix of grid points in the domain specified by x
% FORMAT [X,x] = spm_ndgrid(x)
%
% x{i):   cell array of vectors specifying support or;
% x(i):   vector of bin numbers in the range [-1 1]
%
% x{i):   cell array of vectors specifying support or;
% X:      (n x m) coordinates of n points in m-D space
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ndgrid.m 2032 2008-09-02 18:31:16Z karl $
 
% event-space: domain s
%--------------------------------------------------------------------------
n  = length(x);
if iscell(x)
    s = x;
else
    for i = 1:n
        s{i} = linspace(-1,1,x(i));
    end
end
 
% create X - coordinates of evaluation grid
%---------------------------------------------------------------------------
for i = 1:n
    q     = 1;
    for j = 1:n
        q = spm_kron(s{j}(:).^(i == j),q);
    end
    X(:,i) = q;
end
