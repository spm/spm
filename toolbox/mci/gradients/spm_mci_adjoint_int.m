function [lambda] = spm_mci_adjoint_int (U,P,M,V,djdx,tol)
% Integrate adjoint equation
% FORMAT [lambda] = spm_mci_adjoint_int (U,P,M,V,djdx,tol)
%
% U         Inputs
% P         Parameters
% M         Model structure
% V         states
% djdx      derivative of log likelihood wrt states
% tol       tolerances
%
% lambda    adjoint parameters, at times M.t
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_adjoint_int.m 6548 2015-09-11 12:39:47Z will $

init_t=M.t(1);
final_t=M.T;

tol.abs=1e-3;
tol.rel=1e-5;

% parameters for the intergrator
options = odeset('AbsTol',tol.abs,'RelTol',tol.rel);

% integrate adjoint equation backwards
lam_init=zeros(M.n,1);
[TT, ld] = ode15s(@(t,lam) integrate_adjoints(t,lam,djdx,V,U,P,M),[final_t init_t], lam_init, options);

% interpolate to M.t 
% (Note, interp1q is not accurate enough here
% as TT typically has less resolution than M.t)
lambda = interp1(TT,ld,M.t,'spline');

end

%--------------------------------------------------------------------------
function DlDt = integrate_adjoints(t,lam,djdx,V,U,P,M)

Nx=M.n;
DlDt = zeros(Nx,1);

% Interpolate state and dj/dx to current time point
for d=1:Nx,
    local_djdx(d)=interp1q(M.t,djdx(:,d),t);
    v(d)=interp1q(M.t,V(:,d),t);
end

% Find nearest time point for which we have pre-computed input
if isempty(U)
    ut=[];
else
    [tmp,ind]=min(abs(t-M.t));
    ut=U(:,ind);
end

% Evaluate state Jacobian
if isfield(M,'dfdx')
    Fx = feval(M.dfdx,v(:),ut,P,M);
else
    Fx = spm_diff(M.f,v(:),ut,P,M,1);
end

DlDt=local_djdx-lam'*Fx;
DlDt=DlDt';

end