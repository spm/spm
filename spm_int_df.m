function [y,dydP] = spm_int_df(P,M,U)
% integrates a MIMO nonlinear system using the matrix exponential method,
% and calculates the gradient of the trajectories w.r.t. to system's
% parameters
% FORMAT [x,dydP] = spm_int_df(P,M,U)
% P    - model parameters
% M    - model structure
% U    - input structure or matrix
%
% y    - (v x l) response y = g(x,u,P)
% dydP - partial derivatives of y wrt P
%__________________________________________________________________________
% Integrates the MIMO system described by
%
%    dx/dt = f(x,u,P,M)
%    y     = g(x,u,P,M)
%
% using the matrix exponential (LL) update scheme:
%
%    x(t + dt) = x(t) + Q*f(x,u,P,M)
%
% where the projector Q is defined as:
%
%    (expm(dt*J)-I)*inv(J)
%
% and J is the Jacobian J=df/dx.
%
% This integration scheme also evaluates the partial derivatives (dx/dP) of
% the path x(t) wrt to system's parameters P recursively:
%
%    dx(t+dt)/dP = Q*df/dP + [I+Q*J]*dx(t)/dP + kron(f,I)*dvec(Q)/dP
%
% The interest of the method relies on the speed with which it can
% evaluates dx(t)/dP for all time points. This is based on the analytic
% calculation of the gradients of the flow field w.r.t. the parameters
% (df/dP). This calculation has to be supplied by a function (M.df: see e.g.
% spm_df_erp.m) which outputs both df/dx and df/dP.
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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Jean Daunizeau
% $Id: spm_int_df.m 2353 2008-10-17 11:56:15Z karl $
 
 
% convert U to U.u if necessary
%--------------------------------------------------------------------------
if ~isstruct(U), u.u = U; U = u; end
try, dt = U.dt; catch, dt = 1; end
try, ns = M.ns; catch, ns = length(U.u); end
 
 
% state equation; add [0] states if not specified
%--------------------------------------------------------------------------
try
    f                           = fcnchk(M.f,'x','u','P','M');
catch
    f                           = inline('sparse(0,1)','x','u','P','M');
    M.n                         = 0;
    M.x                         = sparse(0,0);
end
M.nx                            = size(M.x,2);
 
% and output nonlinearity
%--------------------------------------------------------------------------
try
    g                           = fcnchk(M.g,'x','u','P','M');
    % Get expansion point (system's steady state)
    Mog  = rmfield(M,'g');
    [x0] = spm_int_L(P,Mog,U);
    x0   = x0(:,end);
    clear Mog
    % Get derivative of observation function wrt to states at x0
    y0                          = spm_vec(feval(M.g,x0,u,P));
    dgdx0                       = spm_cat(spm_diff(g,x0,u,P,M,1));
catch
    g                           = [];
end
 
% Initial states and inputs
%--------------------------------------------------------------------------
try
    u                           = U.u(1,:);
catch
    u                           = sparse(1,M.m);
end
try
    x                           = M.x;
catch
    x                           = sparse(0,1);
    M.x                         = x;
end
 
% check for evolution function output
%--------------------------------------------------------------------------
try
    dfdx                        = feval(M.df,x,u,P,M,'dfdx');
    dfdP                        = feval(M.df,x,u,P,M,'dfdP');
catch
    error(['flow field derivative function ',M.df,' must return both ',...
        'the Jacobian (df/dx) and derivatives of the flow field wrt ',...
        'system''s parameters (df/dP)!']);
    return
end
 
% get partial derivatives wrt input parameters:
% df/dP_u = (df/du)*(du/dP_u)
dfdP                            = getDfDPu(dfdP,U,1);
% Get projector from Jacobian
[Q,N]                           = getProjector(dfdx,P,M,dt);
 
n                               = numel(x);
np                              = size(dfdP,2);
dxdP                            = zeros(n*(ns+1),np);
if ~isempty(g)
    ny                          = size(y0,1);
    dydP                        = zeros(ny*(ns+1),np);
end
In                              = speye(n);
QJ                              = In+Q*dfdx;
 
 
% Get derivatives of the projector w.r.t. parameters (dQ/dP)
%==========================================================================
% 1- Build an inline function (fQ) that evaluates the projector (Q) from
% the Jacobian (J) that is specified for a given set of parameters
fQ  = @(M,P,x,u,dt,N) getProjector(feval(M.df,x,u,P,M,'dfdx'),P,M,dt,N);
 
% 2- Evaluate numerical derivatives of Q w.r.t. parameters
% !! Get rid of trial effect (B) and input param (R) derivatives, since
% these should be computed outside the ODE integrator (trial effect B is
% computed for each trial, and input time series are computed before this
% function call).
indB                            = spm_getvec(P,'B');
indR                            = spm_getvec(P,'R');
P0                              = rmfield(P,{'B','R'});
dQdp                            = spm_diff(fQ,M,P0,x,u,dt,N,2);
dQdP                            = zeros(n^2,np-length(indR));
for i=1:length(dQdP)
    dQdP(:,i)        = spm_vec(dQdp{i});
end
% 3- Evaluate numerical derivatives of Q w.r.t. first input
[dQdu]                          = spm_diff(fQ,M,P,x,u,dt,N,4);
dQdU = zeros(n^2,length(dQdu));
for i=1:length(dQdu)
    dQdU(:,i)                   = dQdu{i}(:);
end
% 4- Get derivatives w.r.t. input parameters
dQdR                            = getDfDPu(dQdU,U,1);
% 5- Concatenate derivatives wrt ODE params and wrt input params
dQdP                            = [dQdP,dQdR];
 
 
% integrate
%==========================================================================
for i = 1:ns
 
    % input
    %----------------------------------------------------------------------
    try
        u  = U.u(i,:);
    end
 
    % update x(t+dt) = x(t) + dt*f(x(t))
    %----------------------------------------------------------------------
    dxdP(i*n+1:(i+1)*n,:)     = dxdP((i-1)*n+1:i*n,:);
    for j = 1:N
        fx                    = feval(M.f,x,u,P);
        dfdP                  = feval(M.df,x,u,P,M,'dfdP');
        dfdP                  = full(getDfDPu(dfdP,U,i));
        kf                    = kron(fx',In);
        x                     = spm_vec(x) + Q*fx;
        x                     = spm_unvec(x,M.x);
        dxdP(i*n+1:(i+1)*n,:) = Q*dfdP + QJ*dxdP(i*n+1:(i+1)*n,:) + kf*dQdP;
    end
 
 
 
    % output - implement g(x)
    %----------------------------------------------------------------------
    if ~isempty(g)
        y(:,i)                = spm_vec(g(x,u,P,M));
        dydP(i*n+1:(i+1)*n,:) = dgdx0*dxdP(i*n+1:(i+1)*n,:);
    else
        y(:,i)                = spm_vec(x);
    end
 
end
 
% transpose
%--------------------------------------------------------------------------
y                               = real(y');
dxdP(1:n,:)                     = [];
if ~isempty(g)
    dydP(1:n,:)                 = [];
else
    dydP                        = dxdP;
end
 
 
% Sub-functions
%==========================================================================
 
function dfdP = getDfDPu(dfdP,U,t)
% This function gets the partial derivatives w.r.t. the input parameters
% from the partial derivatives w.r.t. inputs time series itself
nu                              = size(U.u,2);
dudP_u                          = permute(U.dudP(t,:,:),[2 3 1]);
dfdu                            = dfdP(:,end-nu+1:end);
dfdP_u                          = dfdu*dudP_u;
dfdP(:,end-nu+1:end)            = [];
dfdP                            = [dfdP,dfdP_u];
 
 
 
function [Q,N] = getProjector(J,P,M,dt,N)
% This function evaluates the projectors required for ODE integration, from
% the Jacobian (J), the parameters (P, which is required for delays) and
% the time discretization dt/N.
nx                              = M.nx;
n                               = size(J,1)./nx;
% propagation delays (intrinsic, extrinsic)
try
    D                           = M.pF.D;
catch % default
    D                           = [2 32];
end
% Delay operator
De                              = D(2).*exp(P.D)/1000;
Di                              = D(1)/1000;
De                              = (1 - speye(n,n)).*De;
Di                              = (1 - speye(nx,nx)).*Di;
De                              = kron(ones(nx,nx),De);
Di                              = kron(Di,speye(n,n));
D                               = Di + De;
iQd                             = speye(nx*n) + D.*J;
Qd                              = pinv(full(iQd));
% Matrix exponential projector
if nargin == 4
    p                           = max(abs(real(eig(full(J)))));
    N                           = ceil(max(1,dt*p*2));
end
iJ                              = pinv(full(J));
Q0                              = expm(dt*Qd*J/N);
Q                               = (Q0 - speye(nx*n,nx*n))*iJ;
% Q = iJ*iQd*(Q0 - speye(9*n,9*n))*Qd;  % this is the LL method
