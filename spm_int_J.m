function [dx] = spm_int_J(J,f,t)
% returns dx(t) = (expm(J*t) - I)*x(0)
% FORMAT [dx] = spm_int_J(J,f,t)
% J      = df/dx
% f      = dx/dt
% dx     = x(t) - x(0)
%--------------------------------------------------------------------------
% Integration of a dynamic system using local linearization.  This scheme
% accommodates nonlinearities in the state equation by using a functional of
% f(x) = dx/dt.  This uses the equality
%
%             expm([0 0]*t) = expm(J*t) - I)*inv(J)*f
%                  [f J]
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id: spm_int_J.m 253 2005-10-13 15:31:34Z guillaume $
 
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
    dx    = spm_int_J(J,f,t);
    
end
 
return
 
% report system is non-dissipative
%--------------------------------------------------------------------------
LE     = max(real(d));
if LE > 0
    warndlg(sprintf('Lyapunov exponent = %.2e',LE))
end




