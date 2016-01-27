function [y,sy,st,x,sx] = spm_mci_sens (P,M,U,csx)
% Integrate dynamics, apply observation model and compute sensitivities
% FORMAT [y,sy,st,x,sx] = spm_mci_sens (P,M,U,csx)
%
% P         Parameters
% M         Model structure
% U         Inputs  [Nin x N]
% csx       Set to 1 to compute state sensitivity
%     
% y         Outputs     [N x Nout]
% sy        Output Sensitivity, dy/dP [N x Nout x Nparams]
% st        Status flag (0 for OK, -1 for problem)
% x         States      [N x Nstates]
% sx        State Sensitivity, dx/dP [N x Nstates x Nparams]
%           ... evaluated at the N time points in M.t
%
% M.f       Flow function dx/dt=f(x,u,P,M)
% M.g       Observation function y=g(x,u,P,M)
%
% This function uses Matlab's ODE suite
%
% B. Sengupta, K. Friston and W. Penny (2014) Efficient Gradient
% Computation for Dynamical Models. Neuroimage,98, 521-527. 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_sens.m 6697 2016-01-27 14:57:28Z spm $

%disp('Warning: spm_mci_sens.m needs fixing for M.V not square !');

y=[];sy=[];x=[];sx=[];
st=0;

try, csx=csx; catch, csx=0; end

% Tolerances for ode15s 
try, tol.rel=M.reltol; catch, tol.rel=1e-2; end
try, tol.abs=M.abstol; catch, tol.abs=1e-4; end

if isstruct(U)
    U=U.u';
end
if isempty(U)
    U=zeros(1,M.N);
end

init_t=0;
final_t=M.T;
init_state=M.x0;

% parameters for the integrator
options = odeset('AbsTol',tol.abs,'RelTol',tol.rel);

% Allocate matrices
x=zeros(M.N,M.n);
y=zeros(M.N,M.l);

try, Np = M.Npflow; catch, Np=length(spm_vec(P)); end
sx=zeros(M.N,M.n,Np);
sy=zeros(M.N,M.l,Np);

% initialise state-sensitivites
Nx=length(M.x0);
init_sens = zeros(Nx*Np,1);

% Use Klopfenstein-Shampine integrator from Matlab's ODE suite
[T,V] = ode15s(@(t,v) integrate_states(t,v,U,P,M),[init_t final_t], [init_state; init_sens], options);

% T     Times
% V     States and sensitivities

% Extract states
xm=V(:,1:Nx);

% Interpolate states to times M.t
x=interp1q(T,xm,M.t);

% Interpolate state sensitivities to times M.t
sxm=interp1q(T,V(:,Nx+1:end),M.t);
sx=reshape(sxm,length(M.t),Nx,Np);

% Compute output and output sensitivity
for n=1:M.N,
    [yout,dydx] = feval (M.g,x(n,:)',U(:,n),P,M);
    y(n,:) = yout'; 
    sy(n,:,:)=dydx*squeeze(sx(n,:,:));
end

end

%--------------------------------------------------------------------------
function DvDt = integrate_states(t,v,U,P,M)
% Integrate states and sensitivities

Nx=length(M.x0);
try, Np = M.Npflow; catch, Np=length(spm_vec(P)); end
Nt=Nx+Nx*Np;  

state_ind=1:Nx;
sens_ind=Nx+1:Nt;

% Find nearest time point for which we have pre-computed input
if isempty(U)
    ut=[];
else
    [tmp,ind]=min(abs(t-M.t));
    ut=U(:,ind);
end

% initialise state matrix
DvDt    = zeros(Nt,1);

% Flow
DvDt(1:Nx)=feval(M.f,v(state_ind),ut,P,M);

if isfield(M,'dfdx')
    Fx = feval(M.dfdx,v(1:Nx),ut,P,M);
else
    Fx = spm_diff(M.f,v(1:Nx),ut,P,M,1);
end

if isfield(M,'dfdp')
    Fp = feval(M.dfdp,v(1:Nx),ut,P,M);
else
    Fp = spm_diff(M.f,v(1:Nx),ut,P,M,3);
end
    
Sx_old = v(sens_ind);
Sx_old = reshape(Sx_old,Nx,Np);

% Gronwall's Theorem
Sx = Fx*Sx_old+Fp;

% Sensitivities
DvDt(sens_ind)=Sx(:);

end



