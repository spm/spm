function [M,stats] = spm_mci_lgv (mcmc,M,U,Y)
% Sampling using Langevin Monte Carlo
% FORMAT [M,stats] = spm_mci_lgv (mcmc,M,U,Y)
%
% mcmc  Sampling parameters
% M     Model Structure
% U     Inputs
% Y     Data
%
% M     Updated model structure
% stats Structure with fields:
%
% .P     Samples, [maxits x M.Np] 
% .E     Negative log joint prob, [maxits x 1]
%
% Uses Simplified Manifold Metropolis Adjusted Langevin 
% Algorithm (Simplified MMALA). 
%
% The manifold matrix captures local curvature but local changes 
% in it are ignored [1,2]. The manifold matrix is more simply 
% interpreted as the posterior covariance under local linear 
% assumptions.
% 
% [1] Calderhead and Girolami. Interface Focus (2011), pp 821-835.
% [2] Girolami and Calderhead. J R Stat Soc B (2011), pp 123-214.
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_lgv.m 6275 2014-12-01 08:41:18Z will $

% Defaults
try verbose=mcmc.verbose; catch verbose=0; end
try maxits=mcmc.maxits; catch maxits=64; end
try plot_int=mcmc.plot_int; catch plot_int=1; end
try update_obs_noise=mcmc.update_obs_noise; catch update_obs_noise=1; end
try update_obs_step=mcmc.update_obs_step; catch update_obs_step=64; end
try h=mcmc.h; catch h=0.5; end 

% Observation noise
Ny=size(Y,2);
noise.c0=Ny;
noise.D0=eye(Ny);

% Compute eigen-parameterisation
M = spm_mci_minit (M);
V  = M.V;

try init=mcmc.init; catch init=M.vpE; end
% Initial param in eigenspace
xinit = M.V'*(init-M.vpE);

% Sample matrix
x = zeros(maxits,M.Np);
x(1,:) = xinit';         

if verbose figure; end

% Tune h by monitoring acceptance rate
tune_h=1;
acc_block=64;
acc_low=0.3;
acc_high=0.7;
total_acc_target=128; % Number of accepted samples to get
acc=zeros(maxits,1);

% Main Loop
i=1; Ce=[];
while (i < maxits) & (sum(acc) < total_acc_target),
    
    if verbose
        if mod(i,plot_int) == 0 & i > 2
            spm_mci_progress (x,E,i);
        end
    end

    if mod(i,acc_block)==0 & tune_h
        % Change step size h ?
        Nacc=sum(acc(i-acc_block+1:i-1));
        prop_acc=Nacc/acc_block;
        if prop_acc < acc_low
            if verbose, disp('Decreasing step size ...'); end
            h=h/2;
        elseif prop_acc > acc_high
            if verbose, disp('Increasing step size ...'); end
            h=h*2;
        end 
    end
    
    if i > 1
        if accepted
            % Only compute gradient and curvature for newly accepted sample
            [j,iCpY,st] = spm_mci_joint_grad (x(i-1,:),M,U,Y);
            if st==-1
                disp('Integration problem in spm_mci_lgv.m');
                keyboard
            end
            
            % Posterior covariance under local linear approximation
            Cp = h*inv(iCpY+M.ipC);
            mu = x(i-1,:)'+0.5*h*Cp*j(:);
        end
        % The above mean and covariance define the proposal density
        % under Langevin dynamics. In MMALA terminology, Cp is the
        % preconditioning matrix. See page 125 of [2].
        
        pos=spm_normrnd(mu,Cp,1);
    else
        pos=x(1,:);
    end
    
    [L,tmp,st] = spm_mci_joint (pos,M,U,Y);
    if st == -1
        disp('Integration problem in spm_mci_lgv.m');
        keyboard
    end
    
    E(i) = -L;
    if i > 1
        dEdit(i-1)=100*(E(i)-E(i-1))/E(i-1);
    end
    [x,E,accepted] = spm_mci_mh_update (x,E,pos,verbose);
    acc(i)=accepted;
    
    % Update observation noise
    if i > update_obs_step & update_obs_noise
        if verbose, disp('Updating observation noise'); end
        try ind=Y.ind; catch ind=1:M.N; end
        if isfield(M,'IS')
            yhat = feval(M.IS,x(i,:)',M,U);
            err=Y-yhat;
        else
            yhat = spm_mci_fwd (x(i,:)',M,U);
            err=Y-yhat(ind,:);
        end
        NT=size(err,1);
        noise.cN=noise.c0+0.5*NT;
        noise.DN=noise.D0+0.5*NT*cov(err);
        Lprec=spm_wishrnd(noise.DN,noise.cN);
        M.Ce=inv(Lprec);
        %M.Ce=iwishrnd(noise.DN,noise.cN);
        Ce(:,:,i-update_obs_step)=M.Ce;
    end
    
    i=i+1;
end

if verbose
    disp(sprintf('Total accepted samples = %d', sum(acc)));
end
% Project parameters back from eigenspace into original space
x=x(1:i-1,:);
nj=size(x,1);
stats.P=M.vpE*ones(1,nj)+V*x';

stats.E=E;
stats.dEdit=dEdit;
stats.acc=acc;
stats.Ce=Ce;