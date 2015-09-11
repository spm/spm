function [P,E,logw,acc] = spm_mci_ais_single_vl (mcmc,M,U,Y,vl)
% Produce a single independent sample from posterior using AIS
% FORMAT [P,E,logw,acc] = spm_mci_ais_single_vl (mcmc,M,U,Y,vl)
%
% mcmc      Sampling settings
% M         Model structure
% U         Input structure
% Y         Data
% vl        Variational Laplace solution
%               .Ep                 Posterior Mean
%               .Cp                 Posterior Covariance
%
% P         [Np x 1] sample
% E         Negative log joint
% logw      Contribution to model evidence
% acc       acc(j) is acceptance rate at temperature j
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_ais_single_vl.m 6548 2015-09-11 12:39:47Z will $

beta=mcmc.beta;
nprop=mcmc.nprop;
prop=mcmc.prop;
scale=mcmc.scale;

J=length(beta);

% Sample from VL posterior (in subspace)
x(:,1) = spm_normrnd (vl.mr,vl.Cr,1);
[L(1),L2] = spm_mci_joint (x(:,1),M,U,Y);

% Sample in original space
xorg = M.V*x(:,1)+M.vpE;
e = xorg - vl.Ep;
log_pvl = vl.const - 0.5*e'*vl.Lambdap*e;
Lsum = (beta(1)-beta(2))*log_pvl+(beta(2)-beta(1))*L(1);

acc(1) = 1;
for j=2:J,
    xs=x(:,j-1); Ls=L(:,j-1);acc(j)=0;
    
    % Generate sample at next temperature using
    % Langevin Monte Carlo
    lgv_mcmc.init=xs;
    lgv_mcmc.maxits=nprop;
    
    [tmp,stats] = spm_mci_lgv_vl (lgv_mcmc,M,U,Y,vl,beta(j));
    
    xs = stats.P(:,end);
    Ls = -stats.E(:,end);
    L2 = stats.L2(end);
    acc (j) = any(stats.acc(2:end));
    
    x(:,j)=xs;
    L(j)=Ls;
    if j < J
        xorg = M.V*x(:,j)+M.vpE;
        e = xorg - vl.Ep;
        log_pvl = vl.const - 0.5*e'*vl.Lambdap*e;
        Lsum = Lsum + (beta(j)-beta(j+1))*log_pvl+(beta(j+1)-beta(j))*L(j);
    end
end
logw=Lsum;
P=x(:,J);
E=-L(J);

if mcmc.rec_traj
    nj=size(x,2);
    traj=spm_vec(M.pE)*ones(1,nj)+V*x;
else
    traj=[];
end