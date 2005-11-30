function [R] = spm_DEM_R(n,s,dt)
% returns the covariance of the temporal derivatives of a Gaussian process
% FORMAT [R] = spm_DEM_R(n,s,dt)
%__________________________________________________________________________
% n    - order of temporal embedding
% s    - temporal smoothness - s.d. of kernel {seconds}
% dt   - time interval {seconds} [default = 1]
%
% R    - (n x n)     E*V*E: covariance of n derivatives
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% defaults
%--------------------------------------------------------------------------
if nargin < 3,        dt = 1; end                 % default sampling time

% no serial dependencies 
%--------------------------------------------------------------------------
if nargin < 2 | ~s, s = 1e-6; end

% temporal correlations (assuming known Gaussian form)
%--------------------------------------------------------------------------
k     = 2*[0:(n - 1)];
x     = sqrt(2)*s/dt;
r(1 + k) = cumprod(1 - k)./x.^k;
R     = [];
for i = 1:n;
    R = [R; r([1:n] + i - 1)];
    r = -r;
end
