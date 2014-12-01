function [MCI] = spm_mci_mfx (MCI)
% Mixed Effects Inference for Dynamical Systems models
% FORMAT [MCI] = spm_mci_mfx (MCI)
%
% MCI               Data structure containing fields:
%
% .M{n}             Model for nth replication (subject)
% .U{n}             Inputs for nth replication
% .Y{n}             Data for nth replication
% .R                Gaussian prior (.pE and .pC) over FFX
% .S                Second level model (optional)  
% .update_ffx       Update fixed effects ? [yes/no] (1/0)
% .update_rfx       Update random effects ? [yes/no] (1/0)
% .update_obs_noise Update observation noise ? [yes/no] (1/0), default=1
% .verbose          Show progress of optimisation 
%
% For Differential Equation models the initial states and flow parameters 
% can be fixed or random effects or take on known values:
%
% .assign.init_par  'fixed', 'random' or 'known'
% .assign.flow_par  'fixed', 'random' or 'known'
%
% .pinit0           Initial values of initial state parameters
%                   [Ninit x 1] for fixed, [Ninit x N] for random
% .pflow0           Initial values of flow parameters
%                   [Nflow x 1] for fixed, [Nflow x N] for random
%
% If the flow_par's are 'known' they should be specified in .pflow0
% If the init_par's are 'known' the M{n}.x0 values are used
%
% The output fields are: 
%
% .sv               [Nv x Nsamples] group fixed effects samples, v
% .sm               [Nw x Nsamples] group random effect means, m
% .sw               [Nw x N x Nsamples] subject random effects, w
% .sv_mean          [Nv x 1] posterior mean over v
% .sm_mean          [Nw x 1] posterior mean over m
% .sw_mean          [Nw x N] posterior mean over w
% .Ce               [Ny x Ny x Nsamples] Obs noise covariance samples
% .postind          Indices for posterior (ie. excluding burn-in)
%
% .noise            Sufficient statistics of conditional distribution
%                   over observation noise precision: .c0,.D0,.cN,.DN
%
% Langevin Monte Carlo implementation described in:
%
% B Sengupta, G Ziegler, G Ridgway, K Friston, J Ashburner and W.Penny.
% Mixed Effects Models of Dynamical Systems, Submitted, 2014.
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_mfx.m 6275 2014-12-01 08:41:18Z will $

try verbose=MCI.verbose; catch verbose=0; end

% Update observation noise ?
try update_obs_noise=MCI.update_obs_noise; catch update_obs_noise=1; end

% Number of samples during Gibbs steps
try mfx_its=MCI.mfx_its; catch mfx_its=256; end
try ffx_its=MCI.ffx_its; catch ffx_its=4; end
try rfx_its=MCI.rfx_its; catch rfx_its=4; end

% Number of RFX/FFX iterations on first MFX iteration
% Note: if fixed_only (see below) ffx1_its is overwritten
% if random_only rfx1_its is overwritten
try ffx1_its=MCI.ffx1_its; catch ffx1_its=64; end
try rfx1_its=MCI.rfx1_its; catch rfx1_its=64; end

M=MCI.M;U=MCI.U;Y=MCI.Y;
Np=length(spm_vec(M{1}.pE));
if isfield(M{1},'x0')
    % For dynamical systems
    d=size(M{1}.x0,1);
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
        Nrand=Nrand+Np;
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
MCI.S=S;

% Observation noise
Ny=size(Y{1}.y,2);
noise.c0=Ny;
noise.D0=eye(Ny);
   
% Prior over fixed effects
try
    fixed=MCI.fixed;
catch
    if strcmp(MCI.assign.init_par,'fixed')
        fixed.pE=MCI.R.pE;
        fixed.pC=MCI.R.pC;
    else
        fixed.pE=[];
        fixed.pC=[];
    end
    if strcmp(MCI.assign.flow_par,'fixed')
        % If flow params are fixed effects their prior mean and cov
        % are set through M{1}.pE, M{1}.pC, not through R.pE, R.pC
        % CHANGE LATER !!
        pE=spm_vec(M{1}.pE);
        pC=M{1}.pC;
        fixed.pE=[fixed.pE;pE];
        fixed.pC=blkdiag(fixed.pC,pC);
    end
    fixed.vpE=spm_vec(fixed.pE); % m_v (vector)
end
MCI.fixed=fixed;
   
% Initialise fixed and random effects
[w_init,v_init,assign] = spm_mci_vw_init (MCI);

% Sampling params
fixed_mcmc.assign=assign;
subj_mcmc.assign=assign;
subj_mcmc.verbose=0;

% Set up sample matrix
Nfixed = length(v_init);
v = zeros(mfx_its,Nfixed);
v(1,:) = v_init';
w = w_init;

% Do we have only RFX or only FFX ?
random_only=1;fixed_only=1;
if strcmp(assign.init_par,'random') | strcmp(assign.flow_par,'random')
    fixed_only=0;
end
if strcmp(assign.init_par,'fixed') | strcmp(assign.flow_par,'fixed')
    random_only=0;
end
if random_only
    disp('Model only has random effects');
    MCI.update_ffx=0;
    rfx1_its=512;
end
if fixed_only 
    disp('Model only has fixed effects');
    MCI.update_rfx=0;
    ffx1_its=512;
    mfx_its=1;
end

% Main Loop
for it=1:mfx_its,
    
    if verbose
        disp(sprintf('MCI iteration %d',it));
    end
        
    if it>1 & Nrand > 0      
        % Update second level params
        S = spm_nwpost (S,w);
        [m,Lambda,C] = spm_nwrnd (S.post,1);
    end
    
    if MCI.update_ffx
        if verbose, disp('Updating fixed effects'); end
        % 2. Update estimate of group fixed effects
        
        % Start sampling where left off
        fixed_mcmc.init=v(it,:)';
        fixed_mcmc.update_obs_noise=update_obs_noise;
        
        if it==1
            fixed_mcmc.verbose=1;
            fixed_mcmc.maxits=ffx1_its;
            if verbose
                disp('Assuming high observation noise precision');
            end
        else
            fixed_mcmc.maxits=ffx_its;
            fixed_mcmc.verbose=0;
        end
        [Psamp,noise,M] = spm_mci_fixed (fixed_mcmc,w,fixed,noise,M,U,Y);
        if fixed_only
            v = Psamp;
        else
            P = Psamp(end,:);
            v(it+1,:) = P;
        end
    else
        v(it+1,:)=v_init';
    end
    
    if MCI.update_rfx
        % Updated prior on random effects
        rfx_prior.pE=m;
        rfx_prior.pC=C;
        
        % 1. Update estimates of subject random effects 
        if verbose, disp('Updating random effects'); end
        if it==1
            subj_mcmc.maxits=rfx1_its;
            subj_mcmc.verbose=0;
            if verbose
                disp('Assuming high observation noise precision');
            end
        else
            subj_mcmc.maxits=rfx_its;
            subj_mcmc.verbose=0;
        end
        
        for n=1:S.N,
            % Start sampling where left off
            subj_mcmc.init=w(:,n);
            
            Psamp = spm_mci_random (subj_mcmc,rfx_prior,v(it+1,:)',M{n},U{n},Y{n});
            w(:,n) = Psamp(:,end);
        end
        
    end
    
    % Update observation noise precision
    if update_obs_noise
        if verbose, disp('Updating observation noise precision'); end
        [noise,M] = spm_mci_obsnoise (w,v(it+1,:)',assign,noise,M,U,Y);
    end
    Ce(:,:,it)=M{1}.Ce;
        
    % Store samples
    sm(:,it)=m;
    sw(:,:,it)=w;
    
end

if fixed_only, total_its=size(v,1); else total_its=mfx_its; end
% Define post burn-in
burn_in=round(0.3*total_its);
post_ind=[burn_in+1:total_its];

% Return samples
MCI.sm=sm;
MCI.sw=sw;
MCI.sv=v';

% Return posterior means
if ~fixed_only
    MCI.sm_mean=mean(sm(:,post_ind),2);
    MCI.sw_mean=mean(sw(:,:,post_ind),3);
end
MCI.sv_mean=mean(v(post_ind,:));
MCI.post_ind=post_ind;

MCI.M=M;
MCI.S=S;
MCI.noise=noise;
MCI.assign=assign;
MCI.Ce=Ce;
