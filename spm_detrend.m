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
% $Id: spm_detrend.m 3696 2010-01-22 14:22:31Z karl $


% defaults
%--------------------------------------------------------------------------
[m n] = size(x);
if ~m || ~n
    y = sparse(m,n);
    return
end
if nargin == 1
    p = 0;
end

% centre columns
%--------------------------------------------------------------------------
if ~p
    for i = 1:n
        y(:,i) = x(:,i) - mean(x(:,i));
    end
    return
end

% polynomial adjustment
%--------------------------------------------------------------------------
G     = [];
for i = 0:p
    d = (1:m).^i;
    G = [G d(:)];
end
y     = x - G*(pinv(full(G))*x);
