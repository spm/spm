function y = spm_cell_swap(x)
% Swap columns for cells in matrix arrays
% FORMAT y = spm_cell_swap(x)
% y{:,i}(:,j) = x{:,j}(:,i);
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% return if empty
%--------------------------------------------------------------------------
if isempty(x)
    y = {};
    return
end

% swap columns for cells
%--------------------------------------------------------------------------
[m,n]  = size(x{1});
[k,l]  = size(x);
y      = cell(k,n);
[y{:}] = deal(zeros(m,l));
for  r = 1:k
    for  j = 1:l
        for  i = 1:n
            y{r,i}(:,j) = x{r,j}(:,i);
        end
    end
end
