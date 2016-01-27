function [dLdp,iCpY,st] = spm_mci_grad_curve (assign,w,v,M,U,Y,fxtype)
% Compute gradient and curvature for MFX model
% FORMAT [dLdp,iCpY,st] = spm_mci_grad_curve (assign,w,v,M,U,Y,fxtype)
%
% assign    fields specify which are random/fixed effects
% w         random effects vector
% v         fixed effects vector
% M,U,Y     structure,inputs,data
% fxtype    'random' or 'fixed'
%
% dLdp      gradient
% iCpY      curvature (Fisher information)
% st        -1 for integration problem
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_grad_curve.m 6697 2016-01-27 14:57:28Z spm $

% Extract init and flow params from rfx or ffx vectors
[p_init,p_flow] = spm_mci_init_flow (assign,w,v,M);

Np=0;
st=0;

% Gradient and curvature of log likelihood
if strcmp(assign.init_par,fxtype)
    [G,sy_init,st] = spm_mci_sens_init (p_init,p_flow,M,U);
    sinit=1;
    Np=Np+size(sy_init,3);
else
    sinit=0;
end

if strcmp(assign.flow_par,fxtype)
    M.x0=p_init; % Initial conditions
    [G,sy_flow,st] = spm_mci_sens (p_flow,M,U);
    sflow=1;
    Np=Np+size(sy_flow,3);
else
    sflow=0;
end

if strcmp(assign.out_par,fxtype)
    M.x0=p_init; % Initial conditions
    % Compute sensitivity to output params
    [G,x,st] = spm_mci_fwd (p_flow,M,U);
    for n=1:M.N,
        [yout,dydx,dydoutp] = feval (M.g,x(n,:)',U(:,n),p_flow,M);
        sy_out(n,:,:)=dydoutp;
    end
    sout=1;
    Np=Np+size(sy_out,3);
else
    sout=0;
end

% Read data points and time indices
try, ind=Y.ind; catch, ind=1:M.N; end
Nt=length(ind);
y=Y.y;

if st==-1, disp('Integration Problem !'); return; end

% Prediction errors
g=G(ind,:);
e=Y.y-g;
       
% Compute gradient and precision
dLdp=zeros(1,Np);
iCpY=zeros(Np,Np);
for t=1:Nt,
    n=ind(t);
    sn=[];
    if sinit
        sn=squeeze(sy_init(n,:,:));
        if M.l==1, sn=sn'; end
    end
    if sflow
        sy=squeeze(sy_flow(n,:,:));
        if M.l==1, sy=sy'; end
        sn=[sn,sy];
    end
    if sout
        sy=squeeze(sy_out(n,:,:));
        if M.l==1, sy=sy'; end
        sn=[sn,sy];
    end
    dLdp=dLdp+e(t,:)*M.iCe*sn;
    iCpY=iCpY+sn'*M.iCe*sn;
end
