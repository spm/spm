function [X] = spm_en(X)
% Euclidean normalization
% FORMAT [X] = spm_en(X);
% X   - matrix
%_______________________________________________________________________
%
% spm_en performs a Euclidean normalization setting column-wise sum of
% squares to unity
%_______________________________________________________________________
% %W% Karl Friston %E%

for i = 1:size(X,2)
	X(:,i) = X(:,i)/sqrt(sum(X(:,i).^2));
end
