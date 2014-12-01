function [dLdp,iCpY] = spm_mci_grad_ind (P,R,M,U,Y,pars)
% Compute gradient and precision wrt selected time points
% FORMAT [dLdp,iCpY] = spm_mci_grad_ind (P,R,M,U,Y,pars)
%
% P         Flow parameters
% R         Initial state parameters
% M         Model structure
% U         Inputs  [Nin x N]
% Y         data; Y.ind denotes selected time points
% pars      'init' or 'flow'
%     
% dLdp      Gradient of log likelihood, dL/dP 
% iCpY      Precision; Inverse Cov of params due to noise on outputs
%
% For pars='init' the function returns gradient and 
% precision wrt initial state parameters. 
%
% For pars='flow' the function returns gradient and 
% precision wrt flow parameters.
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_grad_ind.m 6275 2014-12-01 08:41:18Z will $

% Read data points and time indices
try ind=Y.ind; catch ind=1:M.N; end
Nt=length(ind);
y=Y.y;

switch pars,
    case 'flow',
        M.x0=R; % Initial conditions
        [G,sy,st] = spm_mci_sens (P,M,U);
        dLdp=zeros(1,M.Np);
        iCpY=zeros(M.Np,M.Np);
        
    case 'init',
        [G,sy,st] = spm_mci_sens_init (R,P,M,U);
        dLdp=zeros(1,M.n);
        iCpY=zeros(M.n,M.n);
        
    otherwise
        disp('Unknown parameter type in spm_mci_grad_ind.m');
        return
end

if st==-1, disp('Problem !'); return; end

% Prediction errors
g=G(ind,:);
e=Y.y-g;
       
% Compute gradient and precision
for t=1:Nt,
    n=ind(t);
    sn=squeeze(sy(n,:,:));
    if M.l==1, sn=sn'; end
    dLdp=dLdp+e(t,:)*M.iCe*sn;
    iCpY=iCpY+sn'*M.iCe*sn;
end
