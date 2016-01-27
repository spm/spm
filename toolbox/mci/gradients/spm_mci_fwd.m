function [y,x,st] = spm_mci_fwd (P,M,U)
% Integrate dynamics and apply observation model
% FORMAT [y,x,st] = spm_mci_fwd (P,M,U)
%
% P         Parameters
% M         Model structure
% U         Inputs  [Nin x N]
%     
% y         Outputs     [N x Nout]
% x         States      [N x Nstates]
%           ... evaluated at the N time points in M.t
% st        status flag (0 for OK, -1 for problem)
%
% M.f       Flow function dx/dt=f(x,u,P,M)
% M.g       Observation function y=g(x,u,P,M)
% M.int     Integrator option 
%           eg. 'euler', 'ode15', 'sundials'
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_fwd.m 6697 2016-01-27 14:57:28Z spm $

st=0;
y=[];sy=[];x=[];sx=[];

try, csy=csy; catch, csy=0; end
try, csx=csx; catch, csx=0; end

% Tolerances for ode15s and SUNDIALS
try, reltol=M.reltol; catch, reltol=1e-2; end
try, abstol=M.abstol; catch, abstol=1e-4; end

tDur=[0 M.T];

if isstruct(U)
    U=U.u';
end
if isempty(U)
    U=zeros(1,M.N);
end

switch lower(M.int)
    
    case 'euler',
        x(:,1)=M.x0;
        dt=M.T/M.N;
        for n=2:M.N,
            dxdt=feval(M.f,x(:,n-1),U(:,n-1),P,M);
            x(:,n)=x(:,n-1)+dxdt*dt;
        end
        
    case 'ode15',
        opto = odeset('RelTol',reltol,'AbsTol',abstol);
        sol = ode15s( @(t,x) spm_mci_flow_t(t,x,U,P,M),tDur,M.x0,opto);
        x=deval(sol,M.t);
        
    case 'ode113',
        opto = odeset('RelTol',reltol,'AbsTol',abstol);
        sol = ode113( @(t,x) spm_mci_flow_t(t,x,U,P,M),tDur,M.x0,opto);
        x=deval(sol,M.t);
        
    case 'sundials',
        data.U=U;
        data.P=P;
        data.M=M;
        
        % Refer to the Sundials manual for detailed information on the various integration options
        % e.g. page 8 of sundials/doc/sundialsTB.pdf
        options = CVodeSetOptions('UserData',data,'RelTol',reltol,'AbsTol',abstol, 'LinearSolver','Dense');
        %options = CVodeSetOptions('UserData',data,'RelTol',reltol,'AbsTol',abstol, 'LinearSolver','GMRES');
        CVodeInit(@spm_mci_flow_sun, 'BDF', 'Newton', tDur(1), M.x0, options);
        
            
        [st,time,x] = CVode(M.t,'Normal');
        CVodeFree;
        if st==-1
            disp('Problem in spm_mci_fwd');
            return
        end
        
    otherwise
        disp(sprintf('Error in spm_mci_fwd.m: Unknown integrator %s ',M.int));
        return
end

for n=1:M.N,
    y(:,n) = feval (M.g,x(:,n),U(:,n),P,M);
end

y=y';
x=x';



