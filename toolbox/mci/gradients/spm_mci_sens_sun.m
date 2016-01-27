function [y,sy,st,x,sx] = spm_mci_sens_sun (P,M,U,csx)
% As spm_mci_sens.m but using Sundials
% FORMAT [y,sy,st,x,sx] = spm_mci_sens_sun (P,M,U,csx)
%
% P         Parameters
% M         Model structure
% U         Inputs  [Nin x N]
% csx       Set to 1 to compute state sensitivity
%     
% y         Outputs     [N x Nout]
% sy        Output Sensitivity, dy/dP [N x Nout x Nparams]
% st        Status flag (0 for success, -1 for problem)
% x         States      [N x Nstates]
% sx        State Sensitivity, dx/dP [N x Nstates x Nparams]
%           ... evaluated at the N time points in M.t
%
% M.f       Flow function dx/dt=f(x,u,P,M)
% M.g       Observation function y=g(x,u,P,M)
%
% This function uses the sundials package (CVODE,CVODES,IDA,IDAS)
% from http://computation.llnl.gov/casc/sundials/main.html
%
% B. Sengupta, K. Friston and W. Penny (2014) Efficient Gradient
% Computation for Dynamical Models. Neuroimage,98, 521-527. 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_sens_sun.m 6697 2016-01-27 14:57:28Z spm $

y=[];sy=[];x=[];sx=[];
st=0;

try, csx=csx; catch, csx=0; end

% Tolerances for ode15s and SUNDIALS
try, reltol=M.reltol; catch, reltol=1e-2; end
try, abstol=M.abstol; catch, abstol=1e-4; end

tDur=[0 M.T];
%tDur=[M.t(1) M.T];

if isstruct(U)
    U=U.u';
end
if isempty(U)
    U=zeros(1,M.N);
end

data.U=U;
data.P=P;
data.M=M;

options = CVodeSetOptions('UserData',data,'RelTol',reltol,'AbsTol',abstol, 'LinearSolver','Dense');
CVodeInit(@spm_mci_flow_sun, 'BDF', 'Newton', tDur(1), M.x0, options);

sens_t = zeros(M.n,M.Np);

FSAoptions = CVodeSensSetOptions('method','Simultaneous',...
    'ErrControl', true,...
    'ParamField', 'P',...
    'ParamScales', ones(1,M.Np));

CVodeSensInit(M.Np, [], sens_t, FSAoptions);

time(1,1) = tDur(1);
xm(:,1) = M.x0;
smx(1,:,:) = sens_t;
smy(1,:,:) = zeros(M.l,M.Np);

% When computing output sensitivities, assume dydx=L, ie not a
% function of x. Generalise later
[tmp,L]=feval(M.g,M.x0,U(:,1),P,M);

ti = tDur(1);
it = 1;
while ti < tDur(2)
    it = it+1;
    [status, ti, xt, yS] = CVode(tDur(2),'OneStep');
    if status==-1
        disp('Problem in spm_mci_sens_sun.m');
        st=-1;
        return
    end
    time(it,1) = ti;
    xm(:,it) = xt;
    if csx
        smx(it,:,:) = yS;
    end
    % Edit here if dydx is a function of x
    smy(it,:,:) = L*yS;
end

% Interpolate to M.t
x=interp1(time,xm',M.t,'spline');
x=x';
if csx
    sx=zeros(M.N,M.n,M.Np);
    for j=1:M.Np,
        sx(:,:,j)=interp1q(time,squeeze(smx(:,:,j)),M.t);
    end
end
sy=zeros(M.N,M.l,M.Np);
for j=1:M.Np,
    sy(:,:,j)=interp1q(time,squeeze(smy(:,:,j)),M.t);
end
    
CVodeFree;

for n=1:M.N,
    y(:,n) = feval (M.g,x(:,n),U(:,n),P,M);
end

y=y';
x=x';



