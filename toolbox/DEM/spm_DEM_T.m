function [T] = spm_DEM_T(n,dt)
% returns temporal delay operator
% FORMAT [T] = spm_DEM_T(n,dt)
%__________________________________________________________________________
% n    - order of temporal embedding
% dt   - time interval {seconds} [default = 1]
%
% T    - (n x n)  for generalised state vectors x[t + dt] = T(dt)*x[t]
%
% NB:  T(-dt) = inv(T(dt)), T(-dt)*T(dt) = I and T(i*dT) = T(dt)^i
%==========================================================================
 
% Delay operator (based on Taylor's theorem
%--------------------------------------------------------------------------
T     = zeros(n,n);
for i = 1:n
    t = (dt^(i - 1))/prod(1:(i - 1));
    for j = 1:(n + 1 - i)
        T(j,j + i - 1) = t;
    end
end
