function D = spm_diag_array(X)
% Extract diagonal from 3-D arrays
% FORMAT D = spm_diag_array(X)
%
% X(:,i,i) -> D(:,i);
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2014-2022 Wellcome Centre for Human Neuroimaging


% extract diagnonals
%--------------------------------------------------------------------------
if iscell(X)
    D     = cell(size(X));
    for i = 1:numel(X)
        D{i} = spm_diag_array(X{i});
    end
    return
end
D     = zeros(size(X,1),size(X,2));
for i = 1:size(X,3)
    D(:,i) = X(:,i,i);
end
