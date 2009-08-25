function [x] = spm_cat(x,d)
% converts a cell array into a matrix
% FORMAT [x] = spm_cat(x,d);
% x - cell array
% d - dimension over which to concatenate [default - both]
%__________________________________________________________________________
% Empty array elements are replaced by sparse zero partitions
% and single 0 entries are expanded to conform to the non-empty
% non zero elements.
%
% e.g.:
% > x       = spm_cat({eye(2) []; 0 [1 1; 1 1]})
% > full(x) =
%
%     1     0     0     0
%     0     1     0     0
%     0     0     1     1
%     0     0     1     1
%
% If called with a dimension argument; a cell array is returned
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_cat.m 3333 2009-08-25 16:12:44Z karl $
 
% check x is not already a matrix
%--------------------------------------------------------------------------
if ~iscell(x), return, end
 
% if concatenation over a specific dimension
%--------------------------------------------------------------------------
[n m] = size(x);
if nargin > 1
 
    % concatenate over first dimension
    %----------------------------------------------------------------------
    if d == 1
        y = cell(1,m);
        for i = 1:m
            y{i} = spm_cat(x(:,i));
        end
 
    % concatenate over second
    %----------------------------------------------------------------------
    elseif d == 2
 
        y = cell(n,1);
        for i = 1:n
            y{i} = spm_cat(x(i,:));
        end
 
    % only viable for 2-D arrays
    %----------------------------------------------------------------------
    else
        error('uknown option')
    end
    x      = y;
    return
 
end
 
% find dimensions to fill in empty partitions
%--------------------------------------------------------------------------
for i = 1:n
for j = 1:m
    if iscell(x{i,j})
        x{i,j} = spm_cat(x{i,j});
    end
    [u v]  = size(x{i,j});
    I(i,j) = u;
    J(i,j) = v;
end
end
I     = max(I,[],2);
J     = max(J,[],1);
 
% sparse and empty partitions
%--------------------------------------------------------------------------
[n m] = size(x);
for i = 1:n
for j = 1:m
    if isempty(x{i,j})
        x{i,j} = sparse(I(i),J(j));
    end
end
end
 
% concatenate
%--------------------------------------------------------------------------
for i = 1:n
    y{i,1} = cat(2,x{i,:});
end
try
    x = sparse(cat(1,y{:}));
catch
    x = cat(1,y{:});
end
