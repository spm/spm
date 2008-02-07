function [X] = spm_en(X)
% Euclidean normalization
% FORMAT [X] = spm_en(X);
% X   - matrix
%_______________________________________________________________________
%
% spm_en performs a Euclidean normalization setting column-wise sum of
% squares to unity
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_en.m 1143 2008-02-07 19:33:33Z spm $


for i = 1:size(X,2)
    X(:,i) = X(:,i)/sqrt(sum(X(:,i).^2));
end
