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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_detrend.m 1143 2008-02-07 19:33:33Z spm $


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
