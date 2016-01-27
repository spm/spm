function [S] = spm_mci_random (mcmc,R,v,M,U,Y)
% Random effects estimation
% FORMAT [S] = spm_mci_random (mcmc,R,v,M,U,Y)
%
% mcmc  Sampling parameters
% R     Priors on random effects (R.pE, R.pC)
% v     Fixed effects
% M     Model Structure (single subject)
% U     Inputs (single subject)
% Y     Data (single subject)
%
% S     Samples, [maxits x M.n] 
%
% Uses Langevin Monte Carlo
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_random.m 6697 2016-01-27 14:57:28Z spm $

% Data precision

try, verbose=mcmc.verbose; catch, verbose=0; end
try, maxits=mcmc.maxits; catch, maxits=64; end
try, plot_int=mcmc.plot_int; catch, plot_int=1; end
try, h=mcmc.h; catch, h=0.5; end 

% Assign init/flow/out params as fixed/random effects
try, assign=mcmc.assign; catch, assign=[]; end

% Compute eigen-parameterisation
M = spm_mci_minit (M);

try, init=mcmc.init; catch, init=R.pE; end
% Initial param in eigenspace
xinit = init;

% Read data points and time indices
try, ind=Y.ind; catch, ind=1:M.N; end
Nt=length(ind);
y=Y.y;

% Sample matrix
Nx=size(R.pE,1);
x = zeros(maxits,Nx);
x(1,:) = xinit';          

ipC=pinv(R.pC);
logdet_RCp=spm_logdet(R.pC);

if verbose, figure; end

% Tune h by monitoring acceptance rate
tune_h=1;
acc_block=64;
acc_low=0.3;
acc_high=0.7;
total_acc_target=64; % Number of accepted samples to get
acc=zeros(maxits,1);

i=1;
while (i <= maxits),
    
    if verbose
        if mod(i,plot_int) == 0 && i > 2
            spm_mci_progress (x,E,i);
        end
    end
    
    if mod(i,acc_block)==0 && tune_h
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
    if isfield(M,'IS')
        % Other model types
        [j,iCpY,st,L] = spm_mci_joint_grad (pos,M,U,Y);
    else
        % Differential equation models
        [dLdp,iCpY,st] = spm_mci_grad_curve (assign,pos,v,M,U,Y,'random');
        % Gradient of log prior
        ep = pos-R.pE;
        dlogprior=-ep'*ipC;
        j=dLdp+dlogprior;
        
        [p_init,p_flow] = spm_mci_init_flow (assign,pos,v,M);
        log_like = spm_mci_like_ind (p_flow,p_init,M,U,Y);
        log_prior = -0.5*ep'*ipC*ep + logdet_RCp-0.5*M.Np*log(2*pi);
        L=log_like+log_prior;
    end
    
    if st==-1
        error('Integration problem in spm_mci_random.m');
    end
    prop.L = L;
    
    % Posterior covariance under local linear approximation
    prop.Cp = h*inv(iCpY+ipC);
    prop.mu = pos+0.5*h*prop.Cp*j(:);
    
    prop.iCp = inv(prop.Cp);
    prop.logdetCp = spm_logdet(prop.Cp);
    
    % Accept proposal ?
    [curr,accepted,bayes_fb(i),dL(i)] = spm_mci_mh_update(curr,prop,verbose);
    acc(i)=accepted;
    E(i) = -curr.L;
    if i > 1
        dEdit(i-1)=100*(E(i)-E(i-1))/E(i-1);
    end
    x(i,:)=curr.pos;
    
    i=i+1;
end

S=x(1:i-1,:)';
