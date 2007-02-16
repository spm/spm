function [y] = spm_int_U(P,M,U)
% integrates a MIMO nonlinear system (fast integration for sparse u)
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
% be used.  Otherwise it is evaluated numberically.
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
%--------------------------------------------------------------------------
% %W% Karl Friston %E%

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
for i = 1:M.ns
    
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
    x  = x + Q*fx;
    
    % output - implement g(x)
    %----------------------------------------------------------------------
    if length(g)
        y(:,i) = g(x,u,P,M);
    else
        y(:,i) = x;
    end
    
end

% transpose
%--------------------------------------------------------------------------
y      = real(y');
