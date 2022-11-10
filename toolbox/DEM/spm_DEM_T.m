function [T] = spm_DEM_T(n,dt)
% returns temporal delay operator
% FORMAT [T] = spm_DEM_T(n,dt)
%__________________________________________________________________________
% n    - order of temporal embedding
% dt   - time interval {time steps}
%
% T    - (n x n)  for generalised state vectors x[t + dt] = T(dt)*x[t]
%
% NB:  T(-dt) = inv(T(dt)), T(-dt)*T(dt) = I and T(i*dT) = T(dt)^i
%==========================================================================

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% Delay operator
%--------------------------------------------------------------------------
T = spm_expm(spm_speye(n,n,1)*dt);
    
return

% NOTES

% Delay operator (based on Taylor's theorem)
%--------------------------------------------------------------------------
T     = zeros(n,n);
for i = 1:n
    t = (dt^(i - 1))/prod(1:(i - 1));
    for j = 1:(n + 1 - i)
        T(j,j + i - 1) = t;
    end
end
