function [y] = spm_DEM_embed(Y,n,t,dt)
% temporal embedding into derivatives
% FORMAT [y] = spm_DEM_embed(Y,n,t,dt)
%__________________________________________________________________________
% Y    - (v x N) matrix of v time-series of length N
% n    - order of temporal embedding
% t    - time {secs} at which to evaluate derivatives (starting at t = 1)
% dt   - sampling interval {secs} [default = 1]
%
% y    - {n,1}(v x 1) temporal derivatives   y[:] <- E*Y(t)
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% defaults
%--------------------------------------------------------------------------
if nargin < 4, dt = 1; end

% get dimensions
%--------------------------------------------------------------------------
[q N]  = size(Y);
t      = t/dt;
y      = cell(n,1);
[y{:}] = deal(sparse(q,1));
 
% boundary conditions
%--------------------------------------------------------------------------
k      = [1:n]  + fix(t - (n + 1)/2);
x      = t - min(k) + 1;
i      = k < 1;
k      = k.*~i + i;
i      = k > N;
k      = k.*~i + i*N;


% Inverse embedding operator (T): c.f., a Taylor expansion Y(t) <- T*y[:]
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:n
        T(i,j) = ((i - x)*dt)^(j - 1)/prod(1:(j - 1));
    end
end

% embedding operator: y[:] <- E*Y(t)
%--------------------------------------------------------------------------
E     = inv(T);
 
% embed
%--------------------------------------------------------------------------
for i = 1:n
    y{i} = Y(:,k)*E(i,:)';
end
