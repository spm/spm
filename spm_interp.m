function [x] = spm_interp(x,r)
% 1 or 2-D array interpolation
% FORMAT [x] = spm_interp(x,r)
% x - array
% r - interpolation rate
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_interp.m 5219 2013-01-29 17:07:07Z spm $

% interpolate
%--------------------------------------------------------------------------
[n,m] = size(x);

if n > 1 && m > 1  % matrix
    
    X     = zeros(r*n,m);
    for i = 1:m
        X(:,i) = interpolate(x(:,i),r);
    end
    x     = zeros(r*n,r*m);
    for i = 1:r*n
        x(i,:) = interpolate(X(i,:),r)';
    end

elseif n == 1      % row vector
    x    = interpolate(x',r)';
    
elseif m == 1     % column vector
    x    = interpolate(x,r);

end

% Interpolate using DCT
% -------------------------------------------------------------------------
function [u] = interpolate(y,r)

if r == 1
    u    = y;
else
    y    = y(:);
    n    = size(y,1);
    Dy   = spm_dctmtx(r*n,n);
    Du   = spm_dctmtx(n,n);
    Dy   = Dy*sqrt(r);
    u    = Dy*(Du'*y);
end
