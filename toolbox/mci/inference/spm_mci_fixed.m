function [Psamp,noise,M] = spm_mci_fixed (mcmc,w,fixed,noise,M,U,Y)
% Group fixed effects estimation
% FORMAT [Psamp,noise,M] = spm_mci_fixed (mcmc,w,fixed,noise,M,U,Y)
%
% mcmc      Sampling parameters
% w(:,n)    Random effects for nth subject
% fixed     [fixed.pE, fixed.pC] prior over fixed effects
% noise     noise structure
% M,U,Y     Model, input, data structures
%
% Psamp     Samples, [maxits x M{1}.Np] 
% noise     updated noise model
% M         updated model structures
%
% Uses Langevin Monte Carlo
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_fixed.m 6275 2014-12-01 08:41:18Z will $

try verbose=mcmc.verbose; catch verbose=0; end
try maxits=mcmc.maxits; catch maxits=64; end
try plot_int=mcmc.plot_int; catch plot_int=1; end
try update_obs_step=mcmc.update_obs_step; catch update_obs_step=64; end
try h=mcmc.h; catch h=0.5; end 

assign=mcmc.assign;

% Prior over fixed effects
pE=fixed.pE;
Np=length(pE);
V=eye(Np);
ipC=inv(fixed.pC);
log_prior_t2=spm_logdet(ipC)/2-0.5*Np*log(2*pi);

N=length(M);
for n=1:N
    M{n}.logdet_Ce=spm_logdet(M{n}.Ce);
    M{n}.iCe = inv(M{n}.Ce);
end
    
try init=mcmc.init; catch init=fixed.vpE; end
xinit = init;
x = zeros(maxits,Np);
x(1,:) = xinit';       

if verbose figure; end

% Tune h by monitoring acceptance rate
tune_h=1;
acc_block=32;
acc_low=0.3;
acc_high=0.7;
total_acc_target=64; % Number of accepted samples to get
acc=zeros(maxits,1);

i=1;
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
            % Get gradient and curvature of log like
            dLdp=zeros(1,Np);
            iCpY=zeros(Np,Np);
            Ptmp=x(i-1,:)';
            
            for n=1:N,
                
                if ~isempty(w), wsub=w(:,n); else, wsub=[]; end
                [dLdp_n,iCpY_n] = spm_mci_grad_curve (assign,wsub,Ptmp,M{n},U{n},Y{n},'fixed');
                
                dLdp = dLdp + dLdp_n;
                iCpY = iCpY + iCpY_n;
            end
            
            % Gradient of log prior
            dlogprior = - (x(i-1,:)'-pE)'*ipC;
            
            % Gradient of log joint
            j=dLdp+dlogprior;
            
            % Gaussian proposal
            Cp = h*inv(iCpY+ipC);
            mu = x(i-1,:)'+0.5*h*Cp*j(:);
            
        end
        pos = spm_normrnd(mu,Cp,1);
    else
        pos = x(1,:)';
    end
    Ptmp=pos;
    
    % Log like
    log_like=0;
    for n=1:N,
        if ~isempty(w), wsub=w(:,n); else, wsub=[]; end
        [Pinit,Pflow] = spm_mci_init_flow (assign,wsub,Ptmp,M{n});
        like_n = spm_mci_like_ind (Pflow,Pinit,M{n},U{n},Y{n});
        log_like = log_like + like_n;
    end
    
    % Log prior
    e=pos-pE;
    log_prior = - e'*ipC*e/2 + log_prior_t2;
    
    % Log joint
    L=log_like+log_prior;
    E(i) = -L;
    [x,E,accepted] = spm_mci_mh_update(x,E,pos,verbose);
    acc(i)=accepted;
    
    if mod(i,update_obs_step)==0 & mcmc.update_obs_noise
        if verbose, disp('Updating observation noise'); end
        [noise,M] = spm_mci_obsnoise (w,x(i,:)',assign,noise,M,U,Y);
    end
    i=i+1;
end

Psamp=x(1:i-1,:);
