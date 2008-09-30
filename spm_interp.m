function [x] = spm_interp(x,r)
% 2-D array interpolation
% FORMAT [x] = spm_interp(x,r)
% x - array
% r - interpolation rate
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_interp.m 2271 2008-09-30 21:19:47Z guillaume $

% interpolate
%---------------------------------------------------------------------------
[n m] = size(x);
X     = zeros(r*n,m);
L     = floor((n - 1)/2);
L     = min([L 4]);
for i = 1:m
    X(:,i) = interp(x(:,i),r,L);
end
L     = floor((m - 1)/2);
L     = min([L 4]);
x     = zeros(r*n,r*m);
for i = 1:r*n
    x(i,:) = interp(X(i,:),r,L);
end

