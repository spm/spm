function [y] = spm_cell_swap(x)
% swaps columns for cells in matrix arrays
% FORMAT [y] = spm_cell_swap(x)
% y{:,i}(:,j) = x{:,j}(:,i);
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_cell_swap.m 4310 2011-04-18 16:07:35Z guillaume $

% return if empty
%--------------------------------------------------------------------------
if isempty(x)
    y = {};
    return
end

% swap columns for cells
%--------------------------------------------------------------------------
[m n]  = size(x{1});
[k l]  = size(x);
y      = cell(k,n);
[y{:}] = deal(zeros(m,l));
for  r = 1:k
    for  j = 1:l
        for  i = 1:n
            y{r,i}(:,j) = x{r,j}(:,i);
        end
    end
end



