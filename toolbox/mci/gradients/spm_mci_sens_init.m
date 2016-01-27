function [y,sy,st] = spm_mci_sens_init (R,P,M,U)
% Compute sensitivity to initial state 
% FORMAT [y,sy,st] = spm_mci_sens_init (R,P,M,U)
%
% R         Initial state
% P         Parameters
% M         Model structure
% U         Inputs  [Nin x N]
%     
% y         Outputs     [N x Nout]
% sy        Output Sensitivity, dy/dP [N x Nout x Nparams]
% st        Status flag (0 for OK, -1 for problem)
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
% $Id: spm_mci_sens_init.m 6697 2016-01-27 14:57:28Z spm $

y=[];sy=[];x=[];sx=[];
st=0;

% Tolerances for ode15s 
try, tol.rel=M.reltol; catch, tol.rel=1e-2; end
try, tol.abs=M.abstol; catch, tol.abs=1e-4; end

if isempty(U)
    U=zeros(1,M.N);
end

init_t=0;
final_t=M.T;
init_state=R;

% parameters for the integrator
options = odeset('AbsTol',tol.abs,'RelTol',tol.rel);

% Allocate matrices
x=zeros(M.N,M.n);
y=zeros(M.N,M.l);
sx=zeros(M.N,M.n,M.n);
sy=zeros(M.N,M.l,M.n);

% initialise state-sensitivites
Nx=length(M.x0);
Np=length(P);
S0=eye(M.n);
init_sens = S0(:);

% Use Klopfenstein-Shampine integrator from Matlab's ODE suite
[T,V] = ode15s(@(t,v) int_states_init(t,v,U,P,M),[init_t final_t], [init_state; init_sens], options);

% T     Times
% V     States and sensitivities

% Extract states
xm=V(:,1:Nx);

% Interpolate states to times M.t
x=interp1q(T,xm,M.t);

% Interpolate state sensitivities to times M.t
sxm=interp1q(T,V(:,Nx+1:end),M.t);
%sx=reshape(sxm,length(M.t),Nx,Np);
% Updated Sep 13 2014
sx=reshape(sxm,length(M.t),Nx,Nx);

% When computing output sensitivities, assume dydx=L, ie not a
% function of x. Generalise later
[tmp,L]=feval(M.g,M.x0,U(:,1),P,M);

% Compute output and output sensitivity
for n=1:M.N,
    yout = feval (M.g,x(n,:)',U(:,n),P,M);
    y(n,:) = yout'; 
    sy(n,:,:)=L*squeeze(sx(n,:,:));
end

end

%--------------------------------------------------------------------------
function DvDt = int_states_init (t,v,U,P,M)
% Integrate states and sensitivities to initial conditions

Nx=length(M.x0);
Nt=Nx+Nx*Nx;  

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
DvDt = zeros(Nt,1);

% Flow
DvDt(1:Nx)=feval(M.f,v(state_ind),ut,P,M);

if isfield(M,'dfdx')
    Fx = feval(M.dfdx,v(1:Nx),ut,P,M);
else
    Fx = spm_diff(M.f,v(1:Nx),ut,P,M,1);
end
    
Sx_old = v(sens_ind);
Sx_old = reshape(Sx_old,Nx,Nx);

% Compute change in sensitivity from old - note absence of Fp term
Sx = Fx*Sx_old;

% Sensitivities
DvDt(sens_ind)=Sx(:);

end



