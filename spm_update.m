function [dx] = spm_update(J,f,t)
% returns dx(t) = (expm(J*t) - I)*x(0)
% FORMAT [dx] = spm_update(J,f,t)
% J      = df/dx
% f      = dx/dt
% dx     = x(t) - x(0)
%__________________________________________________________________________


% augment Jacobian and take matrix exponential
%==========================================================================
Jx    = spm_cat({0 []; f J});
dx    = spm_expm(Jx*t);
dx    = dx(2:end,1);

% if system is unstable
%==========================================================================
if norm(dx,1) > 1e6

    % find the eigen-system and remove unstable modes
    %----------------------------------------------------------------------
    [v d] = eig(full(J));
    v     = v(find(real(diag(d))) > 0,:);
    f     = f - v*pinv(v)*f;
    dx    = spm_update(J,f,t);
    
end

return

% report system is non-dissipative
%--------------------------------------------------------------------------
LE     = max(real(d));
if LE > 0
    warndlg(sprintf('Lyapunov exponent = %.2e',LE))
end




