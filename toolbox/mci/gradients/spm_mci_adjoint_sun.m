function [dLdp] = spm_mci_adjoint_sun (Pr,M,U,Y)
% Gradient of log joint from adjoint method (via Sundials)
% FORMAT [dLdp] = spm_mci_adjoint_sun (Pr,M,U,Y)
%
% Pr        Parameters (vectorised and in M.V subspace)
% M         Model structure
% U         Inputs  [Nin x N]
% Y         Data
%     
% dLdp      Gradient [Np x 1]
%
% For M.adjlike=1, dLdp is gradient of log likelihood not log joint
% (useful for debugging).
%
% For M.backint=1 (default), compute the integral underlying dLdp
% *during* backwards integration of adjoint. For M.backint=0, this 
% integral is computed *after* adjoint (useful for debugging).
%
% B. Sengupta, K. Friston and W. Penny (2014) Efficient Gradient
% Computation for Dynamical Models. Neuroimage,98, 521-527. 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_adjoint_sun.m 6697 2016-01-27 14:57:28Z spm $

try, backint=M.backint; catch, backint=1; end
try, adjlike=M.adjlike; catch, adjlike=0; end

% Tolerances 
tol_scale = 1e-3;
reltol = 1e-3;
abstol = 1e-4;

% Parameters in original space
P = M.V*Pr+M.vpE;

%tDur=[0 M.t(end)];
%tDur=[M.t(1) M.t(end)];

t0=0;
%tf=M.T;
%dt=M.t(2)-M.t(1);
%t0=M.t(1)-dt;
tf=M.T;

if isempty(U)
    U=zeros(1,M.N);
end

data.U=U;
data.P=P;
data.M=M;
data.Y=Y;

options = CVodeSetOptions('UserData', data,...
    'RelTol',reltol, ...
    'AbsTol',abstol, ...
    'LinearSolver','Dense', ...
    'JacobianFn',@djacfn);

CVodeInit(@spm_mci_flow_sun, 'BDF', 'Newton', t0, M.x0, options);

CVodeAdjInit(150, 'Hermite');

optionsB = CVodeSetOptions('UserData',data,...
    'MaxNumSteps',50000, ...
    'RelTol',reltol,...
    'AbsTol',abstol,...
    'LinearSolver','Dense',...
    'JacobianFn',@djacBfn);

lambda_init=zeros(M.n,1);
iB = CVodeInitB(@rhsadjoint, 'BDF', ...
    'Newton', tf, lambda_init, optionsB);

if backint
    % Use CVODE to compute the integral underlying dLdp
    % *during* backwards integration of adjoint (backint=1)

    [status,t,y] = CVode(tf,'Normal');
    
    optionsQB = CVodeQuadSetOptions('ErrControl',true,...
        'RelTol',reltol,...
        'AbsTol',abstol);
    
    qB1 = zeros(M.Np,1);
    CVodeQuadInitB(iB, @paramgrad, qB1, optionsQB);
    
    % Backward integration of the adjoint equation
    % (yB is the adjoint vector)
    [status,t,yB,qB] = CVodeB(M.t(2),'Normal');
    %[status,t,yB,qB] = CVodeB(t0,'Normal');
    if status == -1
        error('Adjoint integration failed');
    end
    qB = -qB';
    dt = M.t(2)-M.t(1);
    dLdp = qB/dt;
else
    % Compute integral underlying dLdp *after* backwards
    % integration of adjoint
    
    ntout = 1000;
    dt = (tf-t0)/ntout;
    tt = linspace(t0+dt,tf,ntout-1);
    
    %[status,ttf,x] = CVode(M.t(2:end),'Normal');
    %[status,ttf,x] = CVode(tt,'Normal');
    %[status,ttf,x] = CVode(M.t(2:end),'Normal');
    [status,ttf,x] = CVode(M.t,'Normal');
    interp_x = (interp1q(ttf',x',M.t))';
    
    %tm(1) = tDur(2);
    %t = tDur(2);
    tm(1) = tf;
    t = tf;
    it = 1;
    % Integrate adjoint equation
    while t > M.t(2)
    %while t > t0+2*dt
        it = it+1;
        % The adjoint vector is yB
        %[status,t,yB] = CVodeB(tDur(1),'OneStep');
        [status,t,yB] = CVodeB(M.t(1),'OneStep');
        %[status,t,yB] = CVodeB(t0+2*dt,'OneStep');
        if status == -1
            error('Adjoint integration failed');
        end
        tm(it) = t;
        dldp(it,:) = feval(@paramgrad,t, interp_x, yB, data);
    end
    dldp_int=interp1(tm,dldp,M.t,'spline');
    dLdp=sum(dldp_int);
end
CVodeFree;

if ~adjlike
    dlogpriordp = spm_mci_gprior_deriv (Pr,M);
    dLdp=dLdp+dlogpriordp;
end

end

% -----------------------------------------------------------
function [qBd, flag, new_data] = paramgrad(t, x, yB, data)
% Np quadratures for parametric gradient
% t     time 
% x     state
% yB    adjoint vector, lambda^T
% data  contains P,M,U,Y
%
% yBd   dlambda^T/dt

P = data.P;
M = data.M;
U = data.U;

% If observation function becomes dependent on parameters
% we'll need to update term1
term1 = 0;

% Find nearest time point for which we have pre-computed input
if isempty(U)
    ut=[];
else
    [tmp,ind]=min(abs(t-M.t));
    ut=U(:,ind);
end

% Evaluate parameter Jacobian
if isfield(M,'dfdp')
    dfdp = feval(M.dfdp,x,ut,P,M);
else
    dfdp = spm_diff(M.f,x,ut,P,M,3);
end
term2 = -yB'*dfdp;
qBd = term1 + term2;

flag = 0;
new_data = [];

end

% -----------------------------------------------------------
function [yBd, flag, new_data] = rhsadjoint(t, x, yB, data)
% Adjoint equation
% FORMAT [yBd, flag, new_data] = rhsadjoint(t, x, yB, data)
%
% t     time 
% x     state
% yB    adjoint vector, lambda^T
% data  contains P,M,U,Y
%
% yBd   dlambda^T/dt

Y = data.Y;
M = data.M;
U = data.U;
P = data.P;

% Interpolate data to required time point
for j=1:M.l,
    ydata(j)=interp1q(M.t,Y(:,j),t);
end
[y,L]=feval(M.g,x,U(:,1),P,M);
e=ydata-y';

term1 = (djacfn(t,x,[],data))'*yB; % (f_y)^T \lambda

% When computing output sensitivities, assume dydx=L, ie not a
% function of x. Generalise later
dydx=L;

term2 = e*M.iCe*dydx;
yBd = (-term1 + term2')';

flag = 0;
new_data = [];

end

% -----------------------------------------------------------
function [JB, flag, new_data] = djacBfn(t, y, yB, fyB, data)
% Backward problem Jacobian function

J               = djacfn(t,y,[],data);
JB              = -J';

flag            = 0;
new_data        = [];

end

% -----------------------------------------------------------
function [J, flag, new_data] = djacfn(t, y, fy, data)
% State Jacobian at time t

P=data.P;
M=data.M;
U=data.U;

% Find nearest time point for which we have pre-computed input
if isempty(U)
    ut=[];
else
    [tmp,ind]=min(abs(t-M.t));
    ut=U(:,ind);
end

% Evaluate state Jacobian
if isfield(M,'dfdx')
    J = feval(M.dfdx,y,ut,P,M);
else
    J = spm_diff(M.f,y,ut,P,M,1);
end

flag            = 0;
new_data        = [];

end

