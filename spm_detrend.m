function [y] = spm_detrend(x,p);
% polynomial detrending over columns
% FORMAT [y] = spm_detrend(x,p)
%___________________________________________________________________________
%
% spm_detrend removes linear and nonlinear trends from column-wise data
% matrices
% x   - data matrix
% p   - order of polynomial
% y   - detrended data matrix
%
%__________________________________________________________________________
% %W% %E%

%---------------------------------------------------------------------------
G     = [];
[m n] = size(x);

if ~p
	y = x - ones(m,1)*mean(x);
	return
end
for i = 0:n
	d = [1:m].^p;
	G = [G d(:)];
end

y     = x - G*(G\x);
