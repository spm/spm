function [X] = spm_en(X,p)
% Euclidean normalization
% FORMAT [X] = spm_en(X,[p])
% X   - matrix
% p   - optional polynomial detrend [default: []]
%__________________________________________________________________________
%
% spm_en performs a Euclidean normalization setting the column-wise sum of
% squares to unity (leaving columns of zeros as zeros).
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 1996-2022 Wellcome Centre for Human Neuroimaging


% detrend
%--------------------------------------------------------------------------
if nargin > 1
    X = spm_detrend(X,p);
end

% Euclidean normalization
%--------------------------------------------------------------------------
for i = 1:size(X,2)
    if any(X(:,i))
        X(:,i) = X(:,i)/sqrt(sum(X(:,i).^2));
    end
end
