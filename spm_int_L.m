function [y] = spm_int_L(P,M,U,N)
% Integrate a MIMO nonlinear system using a fixed Jacobian: J(x(0))
% FORMAT [y] = spm_int_L(P,M,U,[N])
% P   - model parameters
% M   - model structure
% U   - input structure or matrix
%
% N   - number of local linear iterations per time step [default: 1]
%
% y   - (v x l)  response y = g(x,u,P)
%__________________________________________________________________________
% Integrates the MIMO system described by
%
%    dx/dt = f(x,u,P,M)
%    y     = g(x,u,P,M)
%
% using the update scheme:
%
%    x(t + dt) = x(t) + U*dx(t)/dt
%
%            U = (expm(dt*J) - I)*inv(J)
%            J = df/dx
%
% at input times.  This integration scheme evaluates the update matrix (U)
% at the expansion point.
%
%--------------------------------------------------------------------------
%
% SPM solvers or integrators
%
% spm_int_ode:  uses ode45 (or ode113) which are one and multi-step solvers
% respectively.  They can be used for any ODEs, where the Jacobian is
% unknown or difficult to compute; however, they may be slow.
%
% spm_int_J: uses an explicit Jacobian-based update scheme that preserves
% nonlinearities in the ODE: dx = (expm(dt*J) - I)*inv(J)*f.  If the
% equations of motion return J = df/dx, it will be used; otherwise it is
% evaluated numerically, using spm_diff at each time point.  This scheme is
% infallible but potentially slow, if the Jacobian is not available (calls
% spm_dx).
%
% spm_int_E: As for spm_int_J but uses the eigensystem of J(x(0)) to eschew
% matrix exponentials and inversion during the integration. It is probably
% the best compromise, if the Jacobian is not available explicitly.
%
% spm_int_B: As for spm_int_J but uses a first-order approximation to J
% based on J(x(t)) = J(x(0)) + dJdx*x(t).
%
% spm_int_L: As for spm_int_B but uses J(x(0)).
%
% spm_int_U: like spm_int_J but only evaluates J when the input changes.
% This can be useful if input changes are sparse (e.g., boxcar functions).
% It is used primarily for integrating EEG models
%
% spm_int:   Fast integrator that uses a bilinear approximation to the
% Jacobian evaluated using spm_bireduce. This routine will also allow for
% sparse sampling of the solution and delays in observing outputs. It is
% used primarily for integrating fMRI models
%__________________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_int_L.m 7143 2017-07-29 18:50:38Z karl $
 
 
% convert U to U.u if necessary
%--------------------------------------------------------------------------
if ~isstruct(U), u.u = U; U = u;   end
try, dt = U.dt;  catch, dt = 1;    end
if nargin < 4; N = 1;              end

% Initial states and inputs
%--------------------------------------------------------------------------
try
    x   = M.x;
catch
    x   = sparse(0,1);
    M.x = x;
end

try
    u   = U.u(1,:);
catch
    u   = sparse(1,M.m);
end

% add [0] states if not specified
%--------------------------------------------------------------------------
try
    f   = spm_funcheck(M.f);
catch
    f   = @(x,u,P,M) sparse(0,1);
    M.n = 0;
    M.x = sparse(0,0);
end
M.f = f;
 
% output nonlinearity, if specified
%--------------------------------------------------------------------------
try
    g = spm_funcheck(M.g);
    if isempty(g)
        g  = @(x,u,P,M) x;
    end
catch
    g = @(x,u,P,M) x;
end
M.g = g;

% dx(t)/dt and Jacobian df/dx (and check for delay operator)
%--------------------------------------------------------------------------
D   = 1;
n   = spm_length(x);
if nargout(f) >= 3
    [fx,dfdx,D] = f(x,u,P,M);
    
elseif nargout(f) == 2
    [fx,dfdx]   = f(x,u,P,M);
    
else
    dfdx = spm_cat(spm_diff(f,x,u,P,M,1)); 
end

% local linear update operator Q = (expm(dt*J) - I)*inv(J)
%--------------------------------------------------------------------------
dfdx  = dfdx - speye(n,n)*exp(-16);
Q     = (spm_expm(dt*D*dfdx/N) - speye(n,n))/dfdx;
 
% integrate
%==========================================================================
v     = spm_vec(x);
for i = 1:size(U.u,1)
    
    % input
    %----------------------------------------------------------------------
    u  = U.u(i,:);
    
    try
        % update dx = (expm(dt*J) - I)*inv(J)*f(x,u)
        %------------------------------------------------------------------
        for j = 1:N
            v = v + Q*f(v,u,P,M);
        end
        
        % output - implement g(x)
        %------------------------------------------------------------------
        y(:,i) = g(v,u,P,M);

    catch
        
        % update dx = (expm(dt*J) - I)*inv(J)*f(x,u)
        %------------------------------------------------------------------
        for j = 1:N
            x = spm_vec(x) + Q*spm_vec(f(x,u,P,M));
            x = spm_unvec(x,M.x);
        end
        
        % output - implement g(x)
        %------------------------------------------------------------------
        y(:,i) = spm_vec(g(x,u,P,M));
 
        
    end
    
end
 
% transpose
%--------------------------------------------------------------------------
y      = real(y');
