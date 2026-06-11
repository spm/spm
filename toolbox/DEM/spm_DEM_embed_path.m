function [y] = spm_DEM_embed_path(Y,n,t,dt)
% Temporal embedding into generalised coordinates of motion
% FORMAT [y] = spm_DEM_embed_path(Y,n,t,dt)
%__________________________________________________________________________
% Y    - (v x N) matrix of v time-series of length N
% n    - order of temporal embedding
% t    - time  {bins} at which to evaluate derivatives (starting at t = 1)
% dt   - sampling interval {secs} [default = 1]
%
% y    - {n,1}(v x 1) temporal derivatives   y[:] <- E*Y(t)
%__________________________________________________________________________
% see also spm_DEM_embed. This version uses slightly different boundary
% conitions

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


% defaults
%--------------------------------------------------------------------------
if nargin < 4, dt = 1; end

% get dimensions
%--------------------------------------------------------------------------
[q,N]  = size(Y);
y      = cell(n,1);
[y{:}] = deal(sparse(q,1));

% return if ~q
%--------------------------------------------------------------------------
if ~q, return, end

% boundary conditions
%--------------------------------------------------------------------------
m  = n;
k  = (1:m)  + fix(t - (m + 1)/2);
x  = t - min(k) + 1;
if max(k) > N
    i = max(k) - N;
    k = k - i;
    x = x + i;
end
k(k < 1) = 1;

% Inverse embedding operator (T): cf, Taylor expansion Y(t) <- T*y[:]
%--------------------------------------------------------------------------
T     = zeros(m,n);
for i = 1:m
    for j = 1:n
        T(i,j) = ((i - x)*dt)^(j - 1)/prod(1:(j - 1));
    end
end

% embedding operator: y[:] <- E*Y(t)
%--------------------------------------------------------------------------
E     = pinv(T);
for i = 1:n
    y{i} = Y(:,k)*E(i,:)';
end
return

