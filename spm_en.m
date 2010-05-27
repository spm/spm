function [X] = spm_en(X,p)
% Euclidean normalization
% FORMAT [X] = spm_en(X,[p]);
% X   - matrix
% p   - optional polynomial detrend [default = []]
%__________________________________________________________________________
%
% spm_en performs a Euclidean normalization setting the column-wise sum of
% squares to unity (leaving columns of zeros as zeros)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_en.m 3901 2010-05-27 16:14:36Z karl $


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
