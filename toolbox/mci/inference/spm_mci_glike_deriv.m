function [dLdp,iCpY,st,L] = spm_mci_glike_deriv (P,M,U,Y)
% Gradient of Gaussian Log-likelihood 
% FORMAT [dLdp,iCpY,st,L] = spm_mci_glike_deriv (P,M,U,Y)
%
% P         Parameters
% M         Model structure
% U         Inputs
% Y         Data
% 
% dLdP      Gradient of Log Likelihood wrt params, [1 x Np]
% iCpY      Curvature (Fisher Information)
% st        status flag (0 for OK, -1 for problem)
% L         Log Likelihood
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_glike_deriv.m 6548 2015-09-11 12:39:47Z will $

%[G,x,S] = spm_mci_fwd (P,M,U,1);
if strcmp(M.int,'sundials')
    [G,S,st] = spm_mci_sens_sun (P,M,U);
else
    [G,S,st] = spm_mci_sens (P,M,U);
end
e = Y-G;

if isfield(M,'Npout') & M.Npout > 0
    % Compute sensitivity to output parameters
    [G,x,st] = spm_mci_fwd (P,M,U);
    sy_out=zeros(M.N,M.l,M.Npout);
    for n=1:M.N,
        [tmp1,tmp2,dydoutp] = feval (M.g,x(n,:)',U(:,n),P,M);
        sy_out(n,:,:)=dydoutp;
    end
    Np=M.Npflow+M.Npout;
else
    Np=size(S,3);
end

dLdp=zeros(1,Np);
for n=1:M.N,
    sn=squeeze(S(n,:,:));
    if M.l==1, sn=sn'; end
    if isfield(M,'Npout') & M.Npout > 0
        sy=squeeze(sy_out(n,:,:));
        if M.l==1, sy=sy'; end
        sn=[sn,sy];
    end
    dLdp=dLdp+e(n,:)*M.iCe*sn;
end

if nargout > 1
    iCpY=zeros(Np,Np);
    for n=1:M.N,
        sn=squeeze(S(n,:,:));
        if M.l==1, sn=sn'; end
        if isfield(M,'Npout') & M.Npout > 0
            sy=squeeze(sy_out(n,:,:));
            if M.l==1, sy=sy'; end
            sn=[sn,sy];
        end
        iCpY=iCpY+sn'*M.iCe*sn;
    end
end

if nargout > 3
    L=spm_mci_glike (P,M,U,Y,G);
end
