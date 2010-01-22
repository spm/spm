function [y] = spm_DEM_embed(Y,n,t,dt,d)
% temporal embedding into derivatives
% FORMAT [y] = spm_DEM_embed(Y,n,t,dt,d)
%__________________________________________________________________________
% Y    - (v x N) matrix of v time-series of length N
% n    - order of temporal embedding
% t    - time  {bins} at which to evaluate derivatives (starting at t = 1)
% dt   - sampling interval {secs} [default = 1]
% d    - delay (bins) for each row of Y
%
% y    - {n,1}(v x 1) temporal derivatives   y[:] <- E*Y(t)
%==========================================================================
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_embed.m 3695 2010-01-22 14:18:14Z karl $

% defaults
%--------------------------------------------------------------------------
if nargin < 4, dt = 1; end
if nargin < 5, d  = 0; end

% get dimensions
%--------------------------------------------------------------------------
[q N]  = size(Y);
y      = cell(n,1);
[y{:}] = deal(sparse(q,1));

% return if ~q
%--------------------------------------------------------------------------
if ~q, return, end

% loop over channels
%--------------------------------------------------------------------------
for p = 1:length(d)

    % boundary conditions
    %----------------------------------------------------------------------
    s      = (t - d(p))/dt;
    k      = [1:n]  + fix(s - (n + 1)/2);
    x      = s - min(k) + 1;
    i      = k < 1;
    k      = k.*~i + i;
    i      = k > N;
    k      = k.*~i + i*N;


    % Inverse embedding operator (T): cf, Taylor expansion Y(t) <- T*y[:]
    %----------------------------------------------------------------------
    for i = 1:n
        for j = 1:n
            T(i,j) = ((i - x)*dt)^(j - 1)/prod(1:(j - 1));
        end
    end

    % embedding operator: y[:] <- E*Y(t)
    %----------------------------------------------------------------------
    E     = inv(T);

    % embed
    %----------------------------------------------------------------------
    if length(d) == q
        for i = 1:n
            y{i}(p,:) = Y(p,k)*E(i,:)';
        end
    else
        for i = 1:n
            y{i}      = Y(:,k)*E(i,:)';
        end
        return
    end
end
