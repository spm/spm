function [X] = spm_en(X)
% Euclidean normalization
% FORMAT [X] = spm_en(X);
% X   - matrix
%_______________________________________________________________________
%
% spm_en performs a Euclidean normalization setting column-wise sum of
% squares to unity
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_en.m 112 2005-05-04 18:20:52Z john $


for i = 1:size(X,2)
	X(:,i) = X(:,i)/sqrt(sum(X(:,i).^2));
end
