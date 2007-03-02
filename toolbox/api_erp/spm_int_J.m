function [y] = spm_int_J(P,M,U)
% integrates a MIMO nonlinear system
% FORMAT [y] = spm_int_J(P,M,U)
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
% at input times.  This integration scheme evaluates the update matrix (Q)
% at each tiem point
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
I   = speye(M.n,M.n);



% integrate
%==========================================================================
for i = 1:ns

    % input
    %----------------------------------------------------------------------
    try
        u  = U.u(i,:);
    end

    % dx(t)/dt and Jacobian df/dx
    %----------------------------------------------------------------------
    try
        [fx dfdx] = f(x,u,P,M);
    catch
        fx        = f(x,u,P,M);
        dfdx      = spm_diff(f,x,u,P,M,1);
    end

    % update dx = (expm(dt*J) - I)*inv(J)*fx
    %----------------------------------------------------------------------
    dx = spm_dx(dfdx,fx,dt);
    x  = x + spm_unvec(dx,x);

    % output - implement g(x)
    %----------------------------------------------------------------------
    if length(g)
        y(:,i) = spm_vec(g(x,u,P,M));
    else
        y(:,i) = spm_vec(x);
    end

end

% transpose
%--------------------------------------------------------------------------
y      = real(y');
