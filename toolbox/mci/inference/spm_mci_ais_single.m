function [P,E,logw,acc,traj] = spm_mci_ais_single (mcmc,M,U,Y)
% Produce a single independent sample using AIS
% FORMAT [P,E,logw,acc,traj] = spm_mci_ais_single (mcmc,M,U,Y)
%
% mcmc      Sampling settings
% M         Model structure
% U         Input structure
% Y         Data
%
% P         [Np x 1] sample
% E         Negative log joint
% logw      Contribution to model evidence
% acc       acc(j) is acceptance rate at temperature j
% traj      traj(p,j) is value of parameter p at temperature j
%           (only set if mcmc.rec_traj=1)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_ais_single.m 6548 2015-09-11 12:39:47Z will $

beta=mcmc.beta;
nprop=mcmc.nprop;
prop=mcmc.prop;
scale=mcmc.scale;

J=length(beta);

if isstruct(M.pC)
    pC = full(diag(spm_vec(M.pC)));
else
    pC = M.pC;
end
% prior cov in subspace
pC=M.V'*pC*M.V;

Np=size(pC,1);

% Sample from prior (in subspace)
x(:,1) = spm_normrnd(zeros(Np,1),pC,1);
[L(1),L2] = spm_mci_joint (x(:,1),M,U,Y);
Lsum = (beta(2)-beta(1))*L2;
acc(1) = 1;
for j=2:J,
    xs=x(:,j-1); Ls=L(:,j-1);acc(j)=0;
    switch prop,
        case 'mh',
            % Generate sample at next temperature by composing
            % Metropolis moves at different scales
            for s=1:nprop,
                dx=spm_normrnd(zeros(Np,1),scale(s)^2*pC,1);
                xcand=xs+dx;
                [Lcand,L2cand] = spm_mci_joint (xcand,M,U,Y,beta(j));
                dL = Lcand-Ls;
                r = exp(dL);
                alpha = min(1,r);
                test_prob = rand(1);
                if alpha > test_prob
                    % Accept
                    xs = xcand;
                    Ls = Lcand;
                    L2 = L2cand;
                    acc (j) = 1;
                end
            end
            
        case 'lmc',
            % Generate sample at next temperature using
            % Langevin Monte Carlo
            M.beta=beta(j);
            lgv_mcmc.init=xs;
            lgv_mcmc.maxits=nprop;
            lgv_mcmc.update_obs_noise=0;
            
            [tmp,stats] = spm_mci_lgv (lgv_mcmc,M,U,Y);
            
            xs = stats.P(:,end);
            Ls = -stats.E(:,end);
            L2 = stats.L2(end);
            acc (j) = any(stats.acc(2:end));
            
        otherwise
            disp('Unknown proposal type in spm_mci_ais_single.m');
    end
    x(:,j)=xs;
    L(j)=Ls;
    if j < J
        Lsum = Lsum + (beta(j+1)-beta(j))*L2;
    end
end
logw=Lsum;
P=x(:,J);
E=-L(J);

if mcmc.rec_traj
    nj=size(x,2);
    traj=spm_vec(M.pE)*ones(1,nj)+M.V*x;
else
    traj=[];
end
