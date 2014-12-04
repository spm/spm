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
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_random.m 6277 2014-12-04 12:16:52Z guillaume $

% Data precision

try, verbose=mcmc.verbose; catch, verbose=0; end
try, maxits=mcmc.maxits; catch, maxits=64; end
try, plot_int=mcmc.plot_int; catch, plot_int=1; end
try, h=mcmc.h; catch, h=0.5; end 

% Assign init/flow params as fixed/random effects
assign=mcmc.assign;

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
while (i < maxits) && (sum(acc) < total_acc_target),
    
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
    
    if i > 1
        if accepted
            if isfield(M,'IS')
                % Other model types
                [dLdp,iCpY] = feval(M.dL,x(i-1,:)',M,U,Y);
            else
                % Differential equation models
                [dLdp,iCpY] = spm_mci_grad_curve (assign,x(i-1,:)',v,M,U,Y,'random');
            end
            
            % Gradient of log prior
            ep = x(i,:)'-R.pE;
            dlogprior=-ep'*ipC;
            
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
    
    % Log like
    if isfield(M,'IS')
        % Other model types
        log_like = feval(M.L,pos,M,U,Y);
    else
        % Differential equation models
        [p_init,p_flow] = spm_mci_init_flow (assign,pos,v,M);
        log_like = spm_mci_like_ind (p_flow,p_init,M,U,Y);
    end
    
    % Log prior
    ep=pos-R.pE;
    %log_prior = -0.5*ep'*ipC*ep + logdet_RCp-0.5*M.n*log(2*pi);
    log_prior = -0.5*ep'*ipC*ep + logdet_RCp-0.5*M.Np*log(2*pi);
    
    % Log joint
    L=log_like+log_prior;
    
    E(i) = -L;
    [x,E,accepted] = spm_mci_mh_update (x,E,pos,verbose);
    acc(i)=accepted;
    
    i=i+1;
end

S=x(1:i-1,:)';
