function [y] = spm_int_U(P,M,U)
% integrates a MIMO nonlinear system (fast integration for sparse inputs)
% FORMAT [y] = spm_int_U(P,M,U)
% P   - model parameters
% M   - model structure
% U   - input structure or matrix
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
% at input times.  This integration scheme is efficient because it
% only evaluates the update matrix (Q) when the inputs change.
% If f returns the Jacobian (i.e. [fx J] = feval(f,M.x,u,P,M) it will
% be used.  Otherwise it is evaluated numerically.
%
% Delay differential equations can be integrated efficiently (but 
% approximately) by absorbing the delay operator into the Jacobian
%
%    dx(t)/dt     = f(x(t - d))
%                 = D(d)f(x(t))
%
%    J(d)         = D(d)df/dx
%
% To invoke this simply delay scheme, f must return D (i.e. [fx J D] =
% feval(f,M.x,u,P,M);
%
% spm_int will also handle static observation models by evaluating
% g(x,u,P,M)
%
% NB: if f returns [] the system is re-set to its initial states.
% This can be useful for implementing boundary conditions.
%
%--------------------------------------------------------------------------
%
% SPM solvers or integrators
%
% spm_int_ode:  uses ode45 (or ode113) which are one and multi-step solvers
% respectively.  They can be used for any ODEs, where the Jacobian is
% unknown or difficult to compute; however, they may be slow.
%
% spm_int_J: uses an explicit Jacobian based update scheme that preserves
% linearities in the ODE: dx = (expm(dt*J) - I)*inv(J)*f.  If the
% equations of motion return J = df/dx, it will be used; otherwise it is
% evaluated numerically, using spm_diff at each time point.  This scheme is
% infallible but potentially slow if the Jacobian is not available (calls
% spm_dx).
%
% spm_int_U: like spm_int_J but only evaluates J when the input changes.
% This can be useful if input changes are sparse (e.g., boxcar functions).
% spm_int_U also has the facility to integrate delay differential equations
% if a delay operator is returned [f J D] = f(x,u,P,M)
%
% spm_int:   Fast integrator that uses a bilinear approximation to the 
% Jacobian evaluated using spm_bireduce. This routine will also allow for
% sparse sampling of the solution and delays in observing outputs
%--------------------------------------------------------------------------
 
% convert U to U.u if necessary
%--------------------------------------------------------------------------
if ~isstruct(U)
    U.u = U;
end
try
    dt = U.dt;
catch
    dt = 1;
end
try
    ns = M.ns;
catch
    ns = length(U.u);
end
 
 
% state equation; add [0] states if not specified
%--------------------------------------------------------------------------
try
    f   = fcnchk(M.f,'x','u','P','M');
catch
    f   = inline('sparse(0,1)','x','u','P','M');
    M.n = 0;
    M.x = sparse(0,0);
end
 
% and output nonlinearity
%--------------------------------------------------------------------------
try
    g   = fcnchk(M.g,'x','u','P','M');
catch
    g   = [];
end
 
% Initial states and inputs
%--------------------------------------------------------------------------
try
    u = U.u(1,:);
catch
    u = sparse(1,M.m);
end
try
    try
        x = feval(M.x0,P,M,U);
    catch
        x = M.x;
    end
catch
    x   = sparse(0,1);
    M.x = x;
end
du  = sparse(1,M.m);
D   = speye(M.n,M.n);
I   = speye(M.n,M.n);
 
% integrate
%--------------------------------------------------------------------------
for i = 1:ns
    
    % input
    %----------------------------------------------------------------------
    try
        du = U.u(i,:) - u;
        u  = U.u(i,:);
    end
    
    % re-compute Jacobian if input changes
    %----------------------------------------------------------------------
    if i == 1 | du*du' > 1e-6
        
        % Jacobian dfdx (evaluated at expansion point)
        %------------------------------------------------------------------
        try
            [fx J D]   = feval(f,M.x,u,P,M);
            J          = D*J;
        catch
            try
                [fx J] = feval(f,M.x,u,P,M);
            catch
                J      = spm_diff(f,M.x,u,P,M,1);
            end
        end
        
        % approximate (expm(dt*J) - I)*inv(J) (avoiding matrix inversion)
        %------------------------------------------------------------------
        T     = I*dt;
        Q     = T;
        for j = 2:256
            T = T*dt*J/j;
            Q = Q + T;
            if norm(T,1) < dt/256, break, end
        end
        Q     = Q*D;
    end
    
    % dx(t)/dt
    %----------------------------------------------------------------------
    fx = f(x,u,P,M);
    
    % reset x(0) if fx = []
    %----------------------------------------------------------------------
    if ~length(fx)
        try
            x = feval(M.x0,P,M,U);
        catch
            x = M.x;
        end
        fx = f(x,u,P,M);
    end
    
    % update dx = (expm(dt*J) - I)*inv(J)*fx = Q*fx;
    %----------------------------------------------------------------------
    x  = x + spm_unvec(Q*fx,x);
    
    % output - implement g(x)
    %----------------------------------------------------------------------
    if length(g)
        y(:,i) = g(x,u,P,M);
    else
        y(:,i) = spm_vec(x);
    end
    
end
 
% transpose
%--------------------------------------------------------------------------
y      = real(y');
