function D = spm_diag_array(X)
% Extracts diagonal from 3-D arrays
% FORMAT D = spm_diag_array(X)
%
% X(:,i,i) -> D(:,i);
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_diag_array.m 5922 2014-03-18 20:10:17Z karl $

% extract diagnonals
%--------------------------------------------------------------------------
if iscell(X)
    D     = cell(size(X));
    for i = 1:length(X)
        D{i} = spm_diag_array(X{i});
    end
    return
end
D     = zeros(size(X,1),size(X,2));
for i = 1:size(X,3);
    D(:,i) = X(:,i,i);
end