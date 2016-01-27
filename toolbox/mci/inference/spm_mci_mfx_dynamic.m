function [MCI] = spm_mci_mfx_dynamic (MCI)
% Mixed Effects Inference for Dynamical Systems
% FORMAT [MCI] = spm_mci_mfx_dynamic (MCI)
%
% MCI               Data structure containing fields:
%
% .M{n}             Model for nth of N replications (e.g. subjects)
% .U{n}             Inputs for nth replication
% .Y{n}             Data for nth replication
% .fixed            Gaussian prior (.pE and .pC) over FFX
% .S                Second level model describing population mean, m, and 
%                   precision, Lambda. The parameters in S.prior
%                   define the sufficient statistics of p(Lambda) (.a and .B)
%                   and p(m|Lambda) (.beta and.m)
% .total_its        Total number of samples per subject
% .rinit            Proportion of samples to collect prior to use of
%                   Empirical (group) prior
% .verbose          Show progress of optimisation
%
% KNOWN, FIXED or RANDOM EFFECTS:
% The initial states, flow and output 
% parameters can be fixed or random effects or take on known values:
%
% .assign.init_par  'fixed', 'random' or 'known'
% .assign.flow_par  'fixed', 'random' or 'known'
% .assign.out_par   'fixed', 'random' or 'known'
%
% .pinit0           Initial values of initial state parameters
%                   [Ninit x 1] for fixed, [Ninit x N] for random
% .pflow0           Initial values of flow parameters
%                   [Nflow x 1] for fixed, [Nflow x N] for random
% .pout0            Initial values of output parameters
%                   [Nout x 1] for fixed, [Nout x N] for random
%
% .update_obs_noise Update observation noise, Gamma ? [yes/no] (1/0), default=1
% .verbose          Show progress of optimisation 
%
% The output fields are: 
%
% POSTERIOR SAMPLES:
% .sv               [Nv x Nsamples] group fixed effects samples, v
% .sm               [Nw x Nsamples] group random effect means, m
% .sw               [Nw x N x Nsamples] subject random effects, w
% .Ce               [Ny x Ny x Nsamples] Obs noise covariance samples
% .postind          Indices for posterior (ie. excluding burn-in)
%
% POSTERIOR MEANS:
% .sv_mean          [Nv x 1] posterior mean over v
% .sm_mean          [Nw x 1] posterior mean over m
% .sw_mean          [Nw x N] posterior mean over w
%
% For Dynamical Systems models:
% .pinit            Estimated initial states
% .pflow            Estimated flow
% .pout             Estimated output params
%
% SUFFICIENT STATISTICS:
% .noise            Parameters of p(Gamma|Y,w,v): .c0,.D0,.cN,.DN
% .S.post           Parameters of p(Lambda|w) (.a and.B) 
%                   and p(m|Lambda,w) (.beta and .m)
%
% W.Penny, M Klein-Flugge and B Sengupta (2015) Mixed Effects Langevin
% Monte Carlo, Submitted, 2015.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mci_mfx_dynamic.m 6697 2016-01-27 14:57:28Z spm $

try, verbose=MCI.verbose; catch, verbose=0; end

% Update observation noise ?
try, update_obs_noise=MCI.update_obs_noise; catch, update_obs_noise=1; end

% Numbers of iterations, with and without group model
try, total_its=MCI.total_its; catch, total_its=1024; end
try, rinit=MCI.rinit; catch, rinit=0.25; end
rfx_its=2;
mfx_its=(1-rinit)*total_its/rfx_its;
rfx1_its=rinit*total_its;

try, ffx_its=MCI.ffx_its; catch, ffx_its=4; end
try, ffx1_its=MCI.ffx1_its; catch, ffx1_its=64; end

M=MCI.M;U=MCI.U;Y=MCI.Y;
Np=length(spm_vec(M{1}.pE));
if isfield(M{1},'x0')
    % For dynamical systems
    d=size(M{1}.x0,1);
end

% Observation noise
Ny=size(Y{1}.y,2);
noise.c0=Ny;
noise.D0=eye(Ny);
   
% Prior over fixed effects
try
    fixed=MCI.fixed;
    fixed.vpE=spm_vec(fixed.pE); % m_v (vector)
catch
    fixed=[];
end
   
% Prior over random effects (second level model)
try
    S=MCI.S;
    Nrand=S.prior.a;
catch
    if strcmp(MCI.assign.init_par,'random')
        Nrand=d;
    else
        Nrand=0;
    end
    if strcmp(MCI.assign.flow_par,'random')
        Nrand=Nrand+M{1}.Npflow;
    end
    if strcmp(MCI.assign.out_par,'random')
        Nrand=Nrand+M{1}.Npout;
    end
    S.prior.P=Nrand;
    S.prior.a=Nrand/2;
    S.prior.B=eye(Nrand);
    S.prior.beta=1;
    S.prior.m=zeros(Nrand,1);
    S.N=length(M);
end
% Prior mean and cov of random effects from second level model
m=S.prior.m;
C=S.prior.B/S.prior.a;
rfx_prior.pE=m;
rfx_prior.pC=C;
MCI.S=S;

% Initialise fixed and random effects
[w_init,v_init,assign,update_ffx,update_rfx] = spm_mci_vw_init (MCI);

% Sampling params
fixed_mcmc.assign=assign;
subj_mcmc.assign=assign;
subj_mcmc.verbose=0;

% Set up sample matrix
Nfixed = length(v_init);
if Nfixed > 0
    %v = zeros(mfx_its,Nfixed);
    for j=1:rfx1_its,
        v(j,:) = v_init';
    end
end
w = w_init;

% Initialisation of each model
subj_mcmc.verbose=verbose;
subj_mcmc.maxits=rfx1_its;
for n=1:S.N,
    subj_mcmc.init=w(:,n);
    if rfx1_its > 0
        Psamp = spm_mci_random (subj_mcmc,rfx_prior,v_init,M{n},U{n},Y{n});
        sw(n,:,:) = Psamp;
        w(:,n) = Psamp(:,end);
    else
        sw(n,:,1) = w(:,n);
    end
end
sw=permute(sw,[2 1 3]);
sm=squeeze(mean(sw,2));
for j=1:rfx1_its,
    Ce(:,:,j)=M{1}.Ce;
end
subj_mcmc.maxits=rfx_its;
subj_mcmc.verbose=0;

% Main Loop
for it=1:mfx_its
    
    if verbose
        disp(sprintf('MCI iteration %d',it));
    end
        
    if it>1 && Nrand > 0   
        % Update second level params
        S = spm_nwpost (S,w);
        [m,Lambda,C] = spm_nwrnd (S.post,1);
    end
    
    
    if update_ffx
        if verbose, disp('Updating fixed effects'); end
        % 2. Update estimate of group fixed effects
        
        % Start sampling where left off
        fixed_mcmc.init=v(it,:)';
        fixed_mcmc.update_obs_noise=update_obs_noise;
        
        if it==1
            fixed_mcmc.verbose=1;
            fixed_mcmc.maxits=ffx1_its;
        else
            fixed_mcmc.maxits=ffx_its;
            fixed_mcmc.verbose=0;
        end
        [Psamp,noise,M] = spm_mci_fixed (fixed_mcmc,w,fixed,noise,M,U,Y);
        P = Psamp(end,:);
        v(rfx1_its+it,:) = P;
        vit=P';
    else
        vit=[];
    end
    
    if update_rfx
        % Updated prior on random effects
        rfx_prior.pE=m;
        rfx_prior.pC=C;
        
        % 1. Update estimates of subject random effects 
        if verbose, disp('Updating random effects'); end
        
        for n=1:S.N,
            if verbose
                disp(sprintf('Subject %d out of %d',n,S.N));
            end
            % Start sampling where left off
            subj_mcmc.init=w(:,n);
            Psamp = spm_mci_random (subj_mcmc,rfx_prior,vit,M{n},U{n},Y{n});
            w(:,n) = Psamp(:,end);
        end
        
    end
    
    % Update observation noise precision
    if update_obs_noise
        if verbose, disp('Updating observation noise precision'); end
        [noise,M] = spm_mci_obsnoise (w,vit,assign,noise,M,U,Y);
    end
    Ce(:,:,rfx1_its+it)=M{1}.Ce;
        
    % Store samples
    sm(:,rfx1_its+it)=m;
    sw(:,:,rfx1_its+it)=w;
    
end

% Define post burn-in
total_its=size(sm,2); 
burn_in=round(0.3*total_its);
post_ind=[burn_in+1:total_its];

% Return samples
MCI.sm=sm;
MCI.sw=sw;
try, MCI.sv=v'; catch, MCI.sv=[]; end

% Return posterior means
if update_rfx
    MCI.sm_mean=mean(sm(:,post_ind),2);
    MCI.sw_mean=mean(sw(:,:,post_ind),3);
end
try, MCI.sv_mean=mean(v(post_ind,:)); catch, MCI.sv_mean=[]; end
MCI.post_ind=post_ind;

% Extract posterior means for pinit, pflow, pout
switch assign.init_par,
    case 'fixed',
        MCI.pinit=MCI.sv_mean(:,assign.v_init)';
    case 'random',
        MCI.pinit_sub=MCI.sw_mean(assign.w_init,:);
        MCI.pinit=MCI.sm_mean(assign.w_init,:);
    otherwise
        try, MCI.pinit=MCI.pinit0; catch, MCI.pinit=[]; end
end
switch assign.flow_par,
    case 'fixed',
        MCI.pflow=MCI.sv_mean(:,assign.v_flow)';
    case 'random',
        MCI.pflow_sub=MCI.sw_mean(assign.w_flow,:);
        MCI.pflow=MCI.sm_mean(assign.w_flow,:);
    otherwise
        try, MCI.pflow=MCI.flow0; catch, MCI.pflow=[]; end
end
switch assign.out_par,
    case 'fixed',
        MCI.pout=MCI.sv_mean(:,assign.v_out)';
    case 'random',
        MCI.pout_sub=MCI.sw_mean(assign.w_out,:);
        MCI.pout=MCI.sm_mean(assign.w_out,:);
    otherwise
        try, MCI.pout=MCI.out0; catch, MCI.pout=[]; end
end

MCI.M=M;
MCI.S=S;
MCI.noise=noise;
MCI.assign=assign;
MCI.Ce=Ce;

