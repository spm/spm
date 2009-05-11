function [X] = spm_en(X)
% Euclidean normalization
% FORMAT [X] = spm_en(X);
% X   - matrix
%__________________________________________________________________________
%
% spm_en performs a Euclidean normalization setting the column-wise sum of
% squares to unity (leaving columns of zeros as zeros)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_en.m 3113 2009-05-11 15:25:13Z karl $

% Euclidean normalization
%--------------------------------------------------------------------------
for i = 1:size(X,2)
    if any(X(:,i))
        X(:,i) = X(:,i)/sqrt(sum(X(:,i).^2));
    end
end
