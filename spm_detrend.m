function y = spm_detrend(x,p)
% Polynomial detrending over columns
% FORMAT y = spm_detrend(x,p)
% x   - data matrix
% p   - order of polynomial [default: 0]
% 
% y   - detrended data matrix
%__________________________________________________________________________
%
% spm_detrend removes linear and nonlinear trends from column-wise data
% matrices.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_detrend.m 5219 2013-01-29 17:07:07Z spm $


% defaults
%--------------------------------------------------------------------------
[m,n] = size(x);
if ~m || ~n
    y = [];
    return
end
if nargin == 1
    p = 0;
end

% centre columns
%--------------------------------------------------------------------------
if ~p
    y = x - ones(m,1)*mean(x);
    return
end

% polynomial adjustment
%--------------------------------------------------------------------------
G     = zeros(m,p+1);
for i = 0:p
    d = (1:m).^i;
    G(:,i+1) = d(:);
end
y     = x - G*(pinv(full(G))*x);
