function [y] = spm_DEM_embed(Y,t,n)
% temporal embedding into derivatives
% FORMAT [y] = spm_DEM_embed(Y,t,n)
%__________________________________________________________________________
% Y    - (v x N) matrix of v time-series of length N
% t    - time (bins) to evaluate derivatives
% n    - order of temporal embedding
%
% y    - {n,1}(v x 1) temporal derivatives   y[:] <- E*Y(t)
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% get dimensions
%--------------------------------------------------------------------------
[m N]  = size(Y);
y      = cell(n,1);
[y{:}] = deal(sparse(m,1));
 
% boundary conditions
%--------------------------------------------------------------------------
k   = [1:n] + t - fix((n + 1)/2);
k   = k((k > 0) & ~(k > N));
n   = length(k);
x   = fix((n + 1)/2);
 
% Inverse embedding operator (T): c.f., a Taylor expansion Y(t) <- T*y[:]
%--------------------------------------------------------------------------
for i = 1:n
    for j = 1:n
        T(i,j) = (i - x)^(j - 1)/prod(1:(j - 1));
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
