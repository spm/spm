function [y] = spm_cell_swap(x)
% swaps columns for cells in matrix arrays
% FORMAT [y] = spm_cell_swap(x)
% y{i}(:,j) = x{j}(:,i);
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_cell_swap.m 272 2005-10-25 20:05:27Z guillaume $

% return if empty
%-----------------------------------------------------------------------
if ~length(x)
    y = {};
    return
end

% swap columns for cells
%-----------------------------------------------------------------------
[m n]  = size(x{1});
l      = length(x);
y      = cell(n,1);
try
    [y{:}] = deal(sparse(m,l));
    for  j = 1:l
        for i = 1:n
            y{i}(:,j) = x{j}(:,i);
        end
    end
end


