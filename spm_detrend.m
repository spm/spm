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
[m n] = size(x);
if nargin == 1, p = 0; end
if ~p
	y = x - ones(m,1)*mean(x);
	return
end
G     = [];
for i = 0:p
	d = [1:m].^i;
	G = [G d(:)];
end

y     = x - G*(pinv(G)*x);
