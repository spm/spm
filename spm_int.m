function [y] = spm_int(P,M,U)
% integrates a MIMO bilinear system dx/dt = f(x,u) = A*x + B*x*u + Cu + D;
% FORMAT [y] = spm_int(P,M,U)
% P   - model parameters
% M   - model structure
% U   - input structure or matrix
%
% y   - (v x l)  response y = g(x,u,P)
%__________________________________________________________________________
% Integrates the bilinear approximation to the MIMO system described by
%
%    dx/dt = f(x,u,P) = A*x + u*B*x + C*u + D
%    y     = g(x,u,P) = L*x;
%
% at v = M.ns is the number of samples [default v = size(U.u,1)]
%
% spm_int will also handle static observation models by evaluating
% g(x,u,P).  It will also handle timing delays if specified in M.delays
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
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id: spm_int.m 1044 2007-12-21 20:36:08Z karl $
 
 
% convert U to U.u if necessary
%--------------------------------------------------------------------------
if ~isstruct(U),
    U.u  = U;
end
try
    U.dt;
catch
    U.dt = 0;
end
 
% get expansion point
%--------------------------------------------------------------------------
x = [1; spm_vec(M.x)];
 
% add [0] states if not specified
%--------------------------------------------------------------------------
if ~isfield(M,'f')
    M.f = inline('sparse(0,1)','x','u','P','M');
    M.n = 0;
    M.x = sparse(0,0);
end
 
% number of times to sample
%--------------------------------------------------------------------------
try
    v = M.ns;
catch
    v = size(U.u,1);
end
 
% output nonlinearity, if specified
%--------------------------------------------------------------------------
try
    g   = fcnchk(M.g,'x','u','P','M');
catch
    g   = inline('x','x','u','P','M');
    M.g = g;
end
 
% Bilinear approximation (1st order)
%--------------------------------------------------------------------------
[M0,M1,L]  = spm_bireduce(M,P);
n          = size(L,2) - 1;                   % n states
m          = size(U.u,2);                     % m inputs
l          = size(L,1);                       % l outputs
u          = size(U.u,1);                     % input times
 
% evaluation time points (when response is sampled or input changes)
%--------------------------------------------------------------------------
if isfield(M, 'delays')
    
    % when delays have been specified transform delays to time bins
    %----------------------------------------------------------------------
    delays = max(1, round(M.delays/U.dt));
    s  = [];
    for j = 1:M.l
        s = [s ceil([0:v-1]*u/v) + delays(j)];
    end
    s  = unique(s);
    s_ind(s) = [1:length(s)]; % index vector from oversampled time to scans
    Nu = length(s);
else
    s  = ceil([1:v]*u/v);     % 'original' output times (last time bin)
    Nu = v;
end
 
t      = [1 (1 + find(any(diff(U.u),2))')];  % input  times
[T s]  = sort([s t]);                        % update (input & output) times
dt     = [U.dt*diff(T) 0];                   % update intervals
 
% Integrate
%--------------------------------------------------------------------------
y      = zeros(l,Nu);
dy     = zeros(l,Nu);
J      = M0;
for  i = 1:length(T)
 
    % input
    %----------------------------------------------------------------------
    u     = U.u(T(i),:);
 
    % change in input - update J
    %----------------------------------------------------------------------
    if s(i) > Nu
 
        J     = M0;
        for j = 1:m
            J = J + u(j)*M1{j};
        end
 
    % output sampled - implement l(x)
    %----------------------------------------------------------------------
    else
        if isfield(M,'g')
            y(:,s(i))  = feval(g,x([1:n] + 1),u,P,M);
        else
            y(:,s(i))  = L*x;
        end
    end
 
    % compute updated states x = expm(J*dt)*x;
    %----------------------------------------------------------------------
    x  = spm_expm(J*dt(i),x);
 
end
 
y      = real(y');
 
% down-sample delays
%----------------------------------------------------------------------
if isfield(M, 'delays')
    
    u    = size(U.u,1);
    tmp  = zeros(v, M.l);
    dtmp = zeros(v, M.l);
    
    % output times for j-th area
    %----------------------------------------------------------------------
    for j = 1:M.l
        s        = ceil([0:v-1]*u/v) + delays(j);
        tmp(:,j) = y(s_ind(s), j);
    end
    y   = tmp;
end
