function [y] = spm_cell_swap(x)
% swaps columns for cells in matrix arrays
% FORMAT [y] = spm_cell_swap(x)
% y{i}(:,j) = x{j}(:,i);
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_cell_swap.m 1143 2008-02-07 19:33:33Z spm $

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


