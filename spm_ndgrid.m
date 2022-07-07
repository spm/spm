function [X,s] = spm_ndgrid(x)
% Return a matrix of grid points in the domain specified by x
% FORMAT [X,x] = spm_ndgrid(x)
%
% x{i):   cell array of vectors specifying support or;
% x(i):   vector of bin numbers in the range [-1 1]
%
% x{i):   cell array of vectors specifying support or;
% X:      (n x m) coordinates of n points in m-D space
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


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
%--------------------------------------------------------------------------
for i = 1:n
    q     = 1;
    for j = 1:n
        q = spm_kron(s{j}(:).^(i == j),q);
    end
    X(:,i) = q;
end
