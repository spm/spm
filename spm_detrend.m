function [y] = spm_detrend(x,p)
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
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_detrend.m 184 2005-05-31 13:23:32Z john $


% defaults
%---------------------------------------------------------------------------
[m n] = size(x);
if ~n
	y = [];
	return
end
if nargin == 1
	p = 0;
end

% centre columns
%---------------------------------------------------------------------------
if ~p
	y = x - ones(m,1)*mean(x);
	return
end

% polynomial adjustment
%---------------------------------------------------------------------------
G     = [];
for i = 0:p
	d = (1:m).^i;
	G = [G d(:)];
end
y     = x - G*(pinv(G)*x);
