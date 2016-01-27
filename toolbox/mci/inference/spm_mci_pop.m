function [P,logev,D,M] = spm_mci_pop (mcmc,M,U,Y)
% Population MCMC with Gaussian priors and proposals
% FORMAT [P,logev,D,M] = spm_mci_pop (mcmc,M,U,Y)
% 
% mcmc      Optimisation parameters  eg.
%
% .J        number of chains 
% .gprob    prob of global move
% .nscale   number of scaling samples
% .ntune    number of tuning samples
% .nsamp    number samples (on avg) to return (per chain)
% .remove_burn_in Remove scale and tune samples.
% .init{j}  [Np x 1] Initial parameter vector for jth chain [optional]
%
% M{i}      Data structure for i=1st or i=2nd model. 1st model
%           is the larger model. Specifying two models is only 
%           necessary if you wish to do model switch integration.
%           Each M structure should contain
%
% .L        Name of log-likelihood function eg. 'spm_dcm_like'
%               must take arguments P,M,U,Y
% .pE       Prior mean
% .pC       Prior covariance
% .lambda1  Observation noise precision
%
% For example, if ith model does not have variable k 
% one can set M{i}.pC(k,k)=1e-4;
%
% U{i}    Input field for ith model (as standard)
% Y       Data field
%
% For each chain the function implements the Adaptive Monte Carlo (AMC)
% algorithm which comprises three phases (i) scaling: proposal cov is 
% optimally scaled prior (same scaling for all params), (ii) tuning: 
% proposal cov is tuned using Robbins-Monro, (iii) sampling: the proposal
% is unchanged. At each stage proposals follow Metropolis-Hastings.
%
% The function returns
%
% P{j}      Posterior samples from jth chain
% logev     approximations to log evidence
%           .pam    Prior Arithmetic Mean 
%           .phm    Posterior Harmonic Mean 
%           .ti     Thermodynamic Integration
%
% For model switch integration logev.ti contains the log Bayes factor
% for model 2 versus model 1
%
% D         Diagnostics
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_pop.m 6697 2016-01-27 14:57:28Z spm $

Nm=length(M);
if Nm>1
    model_switch=1;
else
    model_switch=0;
end
            
try
    J=mcmc.J;
catch
    J=64;
end

try
    gprob=mcmc.gprob;
catch
    %gprob=0.25;
    gprob=0.0;
end

% Precompute some quantities for efficient computation of log prior
if model_switch
    M = spm_mci_switch_prep (M);
else
    M{1} = spm_mci_minit (M{1});
end

V  = M{1}.V;
Np = size(V,1);
Nsub = size(V,2);

if isstruct(M{1}.pC)
    pC = full(diag(spm_vec(M{1}.pC)));
else
    pC = M{1}.pC;
end

try
    anneal=mcmc.anneal;
catch
    anneal='power';
end

% Set annealing schedule
switch anneal
    case 'sigmoid',
        % Good for model switching integration
        x=linspace(-10,10,J);
        beta=1./(1+exp(-x));
    case 'linear',
        beta=linspace(0,1,J);
    case 'log',
        beta=logspace(log10(1/J),0,J-1);
        beta=[0 beta];
    case 'nonlinear',
        eta=0.2;
        j=1:J;
        beta=(eta*j/J)./(1-j/J+eta);
    case 'power',
        % From Calderhead & Girolami
        beta=linspace(0,1,J);
        beta=beta.^5;
    case 'frozen',
        beta=ones(1,J);
end

for j=1:J,
    P{j}.beta=beta(j);
    P{j}.last_update=0;
    P{j}.overall_accept=0;
    P{j}.accept=0;
    P{j}.reject=0;
    P{j}.ar_traj=[];
end

% Initialise proposal covariance
try
    for j=1:J,
        P{j}.C=V'*mcmc.chain{j}.C*V;
    end
catch
    for j=1:J,
        %P{j}.C=V'*M{1}.pC*V;
        P{j}.C=V'*pC*V;
    end
end

% Initialise chains
for j=1:J,
    if isfield(mcmc,'init')
        % User specified
        xinit=mcmc.init{j};
    else
        % Drawn from prior
        if model_switch && isfield(M{1},'Ep')
            xinit=(1-beta(j))*M{1}.Ep+beta(j)*M{2}.Ep;
        else
            xinit=spm_normrnd(M{1}.vpE,M{1}.pC,1);
        end
    end
    P{j}.theta = V'*(spm_vec(xinit)-M{1}.vpE);
    Ep = P{j}.theta;
    if model_switch
        [q,q1,q2] = spm_mci_switch (Ep,M,U,Y,P{j}.beta);
        P{j}.logq(1)=q;
        P{j}.logq1(1)=q1;
        P{j}.logq2(1)=q2;
    else
        [P{j}.logq(1),P{j}.loglike(1)] = spm_mci_joint (Ep,M{1},U{1},Y,P{j}.beta);
    end
end

try
    nscale=mcmc.nscale;
catch
    nscale=500;
end

try
    ntune=mcmc.ntune;
catch
    ntune=500;
end

try
    nsamp=mcmc.nsamp;
catch
    nsamp=1000;
end

try
    remove_burn_in=mcmc.remove_burn_in;
catch
    remove_burn_in=1;
end

% Reset counts for monitoring acceptance rates
for j=1:J,
  P{j}.scale_accepted=0;
  P{j}.scale_rejected=0;
  P{j}.tune_accepted=0;
  P{j}.tune_rejected=0;
  P{j}.scaled=0;
  P{j}.tuned=0;
  P{j}.acc=1;
end

try
    verbose=mcmc.verbose;
catch
    verbose=0;
end

swaps=0;

Ntot=(nscale+ntune+nsamp)*J;
if verbose
    disp('Scaling ...');
end
tic;
dL(1)=0;
for i=2:Ntot,
    
    % Local or Global move
    m=rand(1);
    if m < gprob
        % Global Move
        
        % Now select a pair of chains with adjacent temperatures
        jsel=floor(rand(1)*(J-1))+1;
        
        theta1=P{jsel}.theta(:,end);
        theta2=P{jsel+1}.theta(:,end);
        
        like1=P{jsel}.logq(end);
        like2=P{jsel+1}.logq(end);
        
        Ep1=theta1;
        Ep2=theta2;
        [A,Alike] = spm_mci_joint (Ep2,M,U,Y,P{jsel}.beta);
        [B,Blike] = spm_mci_joint (Ep1,M,U,Y,P{jsel+1}.beta);
        rg=exp(A+B-like1-like2);
        
        % Accept candidate ?
        accept_prob=min(1,rg);

        if rand(1) < accept_prob
            % Accept
            P{jsel}.theta=[P{jsel}.theta theta2];
            P{jsel+1}.theta=[P{jsel+1}.theta theta1];
            P{jsel}.accept=P{jsel}.accept+1;
            P{jsel+1}.accept=P{jsel+1}.accept+1;
            
            P{jsel}.logq(end+1)=A;
            P{jsel+1}.logq(end+1)=B;
            
            P{j}.loglike(end+1)=Alike;
            P{jsel+1}.loglike(end+1)=Blike;
            swaps=swaps+1;
        else
            P{jsel}.theta=[P{jsel}.theta theta1];
            P{jsel+1}.theta=[P{jsel+1}.theta theta2];
            
            P{jsel}.logq(end+1)=P{jsel}.logq(end);
            P{jsel+1}.logq(end+1)=P{jsel+1}.logq(end);
            
            P{jsel}.loglike(end+1)=P{jsel}.loglike(end);
            P{jsel+1}.loglike(end+1)=P{jsel+1}.loglike(end);
        end
        
       
    else
        % Local Move
        
        % Now select a chain
        jsel=floor(rand(1)*J)+1;
        
        if i<nscale*J
            % Step 1 - scale
            P{jsel}.last_update=P{jsel}.last_update+1;
            update_sigma_steps=50;
            if P{jsel}.last_update > update_sigma_steps
                ar=sum(P{jsel}.acc(end-update_sigma_steps+1:end))/update_sigma_steps;
                if ar < 0.2
                    P{jsel}.C=0.5*P{jsel}.C;
                    if verbose
                        disp(sprintf('Acceptance Rate = %1.2f : Decreasing proposal width',ar));
                    end
                elseif ar > 0.4
                    P{jsel}.C=2*P{jsel}.C;
                    if verbose
                        disp(sprintf('Acceptance Rate = %1.2f : Increasing proposal width',ar));
                    end
                else
                    if verbose
                        disp(sprintf('Acceptance Rate = %1.2f',ar));
                    end
                end
                P{jsel}.last_update=0;
            end
            
        elseif i>nscale*J && i < (ntune+nscale)*J
            % Step 2 - tune
            if P{jsel}.scaled==0
                if verbose
                    disp('Tuning ...');
                end
                P{jsel}.scaled=1;
                P{jsel}.adapt_its=1;
                P{jsel}.Ct=P{jsel}.C;
                P{jsel}.mu=mean(P{jsel}.theta,2);
                
                P{jsel}.scale_accepted=P{jsel}.accept;
                P{jsel}.scale_rejected=P{jsel}.reject;
            else
                P{jsel}.adapt_its=P{jsel}.adapt_its+1;
            end
            P{jsel} = spm_mci_update_cov (P{jsel});
            
            if mod(i-nscale,50)==0
                if verbose
                    disp(sprintf('Iteration %d out of %d',i-nscale,ntune));
                end
            end
            
        elseif i>(ntune+nscale)*J
            % Step 3 - sample
            if P{jsel}.tuned==0
                if verbose
                    disp('Sampling ...');
                end
                P{jsel}.tuned=1;
                
                % Eigendecomposition of P{jsel}.Ct 
                sc=2.4/sqrt(Nsub);
               
                if ntune > 0
                    prop_cov=sc^2*P{jsel}.Ct;
                else
                    prop_cov=sc^2*P{jsel}.C;
                end
                [vC,dC]=eig(prop_cov);
                P{jsel}.Cs{1}=diag(dC);
                P{jsel}.Cs{2}=vC;
                
                P{jsel}.tune_accepted=P{jsel}.accept-P{jsel}.scale_accepted;
                P{jsel}.tune_rejected=P{jsel}.reject-P{jsel}.scale_rejected;
            end
            
            if mod(i-nscale-ntune,50)==0
                if verbose
                    disp(sprintf('Iteration %d out of %d',i-nscale-ntune,nsamp));
                end
            end
        end
        
        % Generate candidate
        if i<(nscale+ntune)*J+1
            % Scaling and tuning phases
            sc=2.4/sqrt(Nsub);
            prop_cov=sc^2*P{jsel}.C;
            xcand=P{jsel}.theta(:,end)+sqrt(diag(prop_cov)).*randn(Nsub,1);
        else
            % Sampling phase
            xcand=P{jsel}.theta(:,end)+spm_normrnd(zeros(Nsub,1),P{jsel}.Cs,1);
        end
        
        % Candidates are already in reduced space
        Epcand = xcand;
        if model_switch
            [cand_like,logq1,logq2] = spm_mci_switch (Epcand,M,U,Y,P{jsel}.beta);
        else
            [cand_like,loglike] = spm_mci_joint (Epcand,M{1},U{1},Y,P{jsel}.beta);
        end
      
        % Accept candidate ?
        dL(i)=cand_like-P{jsel}.logq(end);
        accept_prob=min(1,exp(dL(i)));
        
        if rand(1) < accept_prob
            % Accept
            P{jsel}.theta=[P{jsel}.theta xcand];
            P{jsel}.logq(end+1)=cand_like;
            P{jsel}.accept=P{jsel}.accept+1;
            P{jsel}.overall_accept=P{jsel}.overall_accept+1;
            
            if model_switch
                P{jsel}.logq1(end+1)=logq1;
                P{jsel}.logq2(end+1)=logq2;
            else
                P{jsel}.loglike(end+1)=loglike;
            end
            
            P{jsel}.acc=[P{jsel}.acc;1];
            
        else
            % Reject
            P{jsel}.theta(:,end+1)=P{jsel}.theta(:,end);
            P{jsel}.logq(end+1)=P{jsel}.logq(end);
            P{jsel}.reject=P{jsel}.reject+1;
            
            if model_switch
                P{jsel}.logq1(end+1)=P{jsel}.logq1(end);
                P{jsel}.logq2(end+1)=P{jsel}.logq2(end);
            else
                P{jsel}.loglike(end+1)=P{jsel}.loglike(end);
            end
            
            P{jsel}.acc=[P{jsel}.acc;0];
        end
        
    end
end


burn_in=nscale+ntune;
if remove_burn_in
    % Remove burn-in samples from each chain
    for j=1:J,
        P{j}.theta=P{j}.theta(:,burn_in+1:end);
        P{j}.logq=P{j}.logq(burn_in+1:end);
        if model_switch
            P{j}.logq1=P{j}.logq1(burn_in+1:end);
            P{j}.logq2=P{j}.logq2(burn_in+1:end);
        else
            P{j}.loglike=P{j}.loglike(burn_in+1:end);
        end
    end
end

% Project parameters from eigenspace of prior back into original space
for j=1:J,
    nj=size(P{j}.theta,2);
    P{j}.theta=spm_vec(M{1}.pE)*ones(1,nj)+V*P{j}.theta;
end
    
for j=1:J,
    P{j}.dL=dL;
    P{j}.ar=P{j}.overall_accept/(burn_in+nsamp);
    D.ar(j)=P{j}.ar;
    
    P{j}.ar_scale=P{j}.scale_accepted/(P{j}.scale_accepted+P{j}.scale_rejected);
    P{j}.ar_tune=P{j}.tune_accepted/(P{j}.tune_accepted+P{j}.tune_rejected);            
    
    sample_accepted=P{j}.overall_accept-P{j}.scale_accepted-P{j}.tune_accepted;
    sample_rejected=P{j}.reject-P{j}.scale_rejected-P{j}.tune_rejected;
    P{j}.ar_sample=sample_accepted/(sample_accepted+sample_rejected);
    beta(j)=P{j}.beta;
    
    if model_switch
        Uhat(j)=mean(P{j}.logq2-P{j}.logq1);
    else
        Uhat(j)=mean(P{j}.loglike);
    end
end

% Compute evidence using trapezoidal rule:
% Lartillot and Phillipe (2006) Computing Bayes Factors Using Thermodynamic
% Integration, Systematic Biology, 2006, 55(2):195-207.
if J>1
    logev.ti=trapz(beta,Uhat);
else
    logev.ti=Uhat;
end

D.swaps=swaps;
D.Uhat=Uhat;
D.beta=beta;

% For frozen annealing (ie plain MH) compute overall PHM
if strcmp(anneal,'frozen')
    L=[];Lq=[];
    for j=1:J,
        L=[L, P{j}.loglike];
        Lq=[Lq, P{j}.logq];
    end
    logev.phm = spm_mci_phm (L);
else
    logev.phm = [];
end

D.els=toc;


