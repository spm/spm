function D = spm_diag_array(X)
% Extracts diagonal from 3-D arrays
% FORMAT D = spm_diag_array(X)
%
% X(:,i,i) -> D(:,i);
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_diag_array.m 5956 2014-04-16 14:34:25Z guillaume $

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
