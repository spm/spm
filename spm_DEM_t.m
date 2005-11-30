function [pt] = spm_DEM_t(n,r,s,dt)
% FORMAT [pt] = spm_DEM_t(n,r,s,dt)
% FORMAT [pt] = spm_DEM_t(M)
%__________________________________________________________________________
% n    - embedding order 
% d    - restriction order
% s    - temporal smoothness - s.d. of kernel {seconds}
% dt   - time interval {seconds}
%
% pt  - prediction interval = max F(pt) = min KL{P(0,n) P(pt,d)}
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% get parameters from model structure
%--------------------------------------------------------------------------
if isstruct(n)
    dt   = n(1).E.dt;                        % time step
    r    = n(1).E.r;                         % restriction order
    s    = n(1).E.s;                         % smoothness - s.d. of kernel
    n    = n(1).E.n;                         % order of embedding (n >= d)
end

% if there is no temporal embedding (n = 1), pt = 1
%--------------------------------------------------------------------------
if n == 1
    pt = 1;
    return
end

if nargout
    pt = fminbnd(@(pt) spm_DEM_t_fmin(n,r,s,dt,pt),0,8);
    return
end

% line search over F - D{P(pt)||P(0)) for graphics
%==========================================================================

% prediction times
%--------------------------------------------------------------------------
N     = 128;
PT    = s*8*[1:N]/N;
for i = 1:N
    P    = spm_DEM_P(n,r,s,dt,PT(i));
    L    = svd(full(P));
    F(i) = log(prod(L(1:r)))/2;
end

% estimated pt = max(F(pt))
%--------------------------------------------------------------------------
[j i] = max(real(F));
pt    = PT(i);

% graphics
%==========================================================================
plot(PT,F,[pt pt],[min(F) max(F)],':')
ylabel('- Free energy')
xlabel('Interval {secs}')
axis square

