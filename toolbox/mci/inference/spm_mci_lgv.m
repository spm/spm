function [M,stats] = spm_mci_lgv (mcmc,M,U,Y)
% Sampling using Langevin Monte Carlo
% FORMAT [M,stats] = spm_mci_lgv (mcmc,M,U,Y)
%
% mcmc  Sampling parameters
%       .verbose            display progress
%       .maxits             maximum number of total samples 
%       .init               initial sample values (start of chain)
%       .update_obs_noise   estimate observation noise
%       .update_obs_step    update obs noise after this number of samples
%       .restart            restart chain from init 
%       .h                  step size
%       .adapt_h            adapt h based on acceptance rate
%
% M     Model Structure
%       .dL                 Gradients and curvatures are computed using 
%                           this user-specified function. If this is absent
%                           they will be computed using (i) the forward
%                           sensitivity method for dynamical models 
%                           (ie. if M.f exists) or (ii) finite differences
%                           otherwise
%                           
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
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_lgv.m 6697 2016-01-27 14:57:28Z spm $

% Defaults
try, verbose=mcmc.verbose; catch, verbose=0; end
try, maxits=mcmc.maxits; catch, maxits=64; end
try, restart=mcmc.restart; catch, restart=0; end
try, plot_int=mcmc.plot_int; catch, plot_int=1; end
try, update_obs_noise=mcmc.update_obs_noise; catch, update_obs_noise=1; end
try, update_obs_step=mcmc.update_obs_step; catch, update_obs_step=maxits/2; end
try, h=mcmc.h; catch, h=0.5; end 
try, adapt_h=mcmc.adapt_h; catch, adapt_h=0; end
try, init=mcmc.init; catch, init=spm_vec(M.pE); end

% Observation noise
Ny=size(Y,2);
noise.c0=Ny;
noise.D0=eye(Ny);

% Compute eigen-parameterisation
M = spm_mci_minit (M);
V  = M.V;

% Initial param in eigenspace
xinit = M.V'*(init-M.vpE);

% Sample matrix
x = zeros(maxits,M.Np);
x(1,:) = xinit';         

if verbose figure; end

% Parameters for adapting h by monitoring acceptance rate
acc_block=64;
acc_low=0.3;
acc_high=0.7;
    
acc=zeros(maxits,1);

if ~isfield(M,'Ce')
    % Don't update observation noise if we don't have a 
    % Gaussian noise model
    update_obs_noise=0;
end

% Main Loop
i=1; Ce=[];
while (i <= maxits),
    
    if verbose
        if mod(i,plot_int) == 0 && i > 2
            spm_mci_progress (x,E,i);
        end
    end

    if mod(i,acc_block)==0 && adapt_h
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
    
    % Proposal (first proposal always accepted)
    if i==1
        pos=x(1,:)';
        curr=[];
    else
        pos=spm_normrnd(curr.mu,curr.Cp,1);
    end
        
    % Quantities re proposal
    prop.pos=pos;
    [j,iCpY,st,prop.L,prop.L2] = spm_mci_joint_grad (pos,M,U,Y);
    if st==-1
        error('Integration problem in spm_mci_lgv.m');
    end
    % Posterior covariance under local linear approximation
    prop.Cp = h*inv(iCpY+M.ipC);
    prop.mu = pos+0.5*h*prop.Cp*j(:);
    
    prop.iCp = inv(prop.Cp);
    prop.logdetCp = spm_logdet(prop.Cp);
    
    % Accept proposal ?
    [curr,accepted,bayes_fb(i),dL(i)] = spm_mci_mh_update(curr,prop,verbose);
    acc(i)=accepted;
    E(i) = -curr.L;
    L2(i) = curr.L2;
    if i > 1
        dEdit(i-1)=100*(E(i)-E(i-1))/E(i-1);
    end
    x(i,:)=curr.pos;
    
    % Update observation noise
    if (i > update_obs_step && update_obs_noise) || (restart && update_obs_noise)
        if verbose, disp('Updating observation noise'); end
        try, ind=Y.ind; catch, ind=1:M.N; end
        % Parameters in original space
        xorg = M.V*x(i,:)'+M.vpE;
        if isfield(M,'IS')
            %yhat = feval(M.IS,x(i,:)',M,U);
            yhat = feval(M.IS,xorg,M,U);
            err=Y-yhat;
        else
            %yhat = spm_mci_fwd (x(i,:)',M,U);
            yhat = spm_mci_fwd (xorg,M,U);
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
stats.bayes_fb=bayes_fb;
stats.dL=dL;
stats.L2=L2;
