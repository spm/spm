function [R,V,D] = spm_DEM_R(n,s,dt)
% returns the precision of the temporal derivatives of a Gaussian process
% FORMAT [R,V,D] = spm_DEM_R(n,s,[dt])
%__________________________________________________________________________
% n    - truncation order
% s    - temporal smoothness - s.d. of kernel {seconds}
% dt   - bin interval (seconds)
%
%                         e[:] <- E*e(0)
%                         e(0) -> D*e[:]
%                 <e[:]*e[:]'> = R
%                              = <E*e(0)*e(0)'*E'>
%                              = E*V*E'
%
% R    - (n x n)     E*V*E: precision of n derivatives
% V    - (n x n)     V:    covariance of n derivatives
% D    - (n x n)     D:    projector:   Y(t) <- D*y[:]
%==========================================================================

% if no serial dependencies 
%--------------------------------------------------------------------------
if ~n,              R = sparse(0,0);
                    V = R;return, end
if nargin < 2 | ~s, s = exp(-16); end
if nargin < 3,     dt = 1;        end


% temporal correlations (assuming known Gaussian form) - V
%--------------------------------------------------------------------------
k     = 2*[0:(n - 1)];
x     = sqrt(2)*s;
r(1 + k) = cumprod(1 - k)./x.^k;
V     = [];
for i = 1:n;
    V = [V; r([1:n] + i - 1)];
    r = -r;
end

% precision - R
%--------------------------------------------------------------------------
R     = inv(V);

% Inverse embedding operator (D): c.f., a Taylor expansion Y(t) <- D*y[:]
%--------------------------------------------------------------------------
x     = fix((n + 1)/2);
for i = 1:n
    for j = 1:n
        D(i,j) = ((i - x)*dt)^(j - 1)/prod(1:(j - 1));
    end
end

if nargout, return, end

% graphics
%==========================================================================

% embedding operator: y[:] <- E*Y(t): Graphics
%--------------------------------------------------------------------------
t     = ([1:n] - x)*dt;
subplot(2,2,1)
imagesc(V)
axis square
title({'covariance';'derivatives'})

subplot(2,2,2)
imagesc(t,t,D*V*D')
axis square
title({'covariance';'time (secs)'})

subplot(2,2,3)
imagesc(R)
% plot(R)
axis square
title({'precision';'derivatives'})

subplot(2,2,4)
imagesc(t,t,inv(D*V*D'))
% plot(inv(D*V*D'))
axis square
title({'precision';'time (secs)'})

return

% NB: temporal correlations (assuming stationary Gaussian form)
%--------------------------------------------------------------------------
t     = ([1:n] - 1)*dt;
v     = 2*(s^2);
V     = exp(-t.^2/(2*v));
V     = toeplitz(V);


