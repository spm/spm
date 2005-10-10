function [y] = spm_cell_swap(x)
% swaps columns for cells in matrix arrays
% FORMAT [y] = spm_cell_swap(x)
% y{i}(:,j) = x{j}(:,i);
%__________________________________________________________________________
% %W% Karl Friston %E%

% return if empty
%--------------------------------------------------------------------------
if ~length(x)
    y = {};
    return
end

% swap columns for cells
%--------------------------------------------------------------------------
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


