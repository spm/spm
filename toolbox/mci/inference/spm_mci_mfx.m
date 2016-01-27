function [MCI] = spm_mci_mfx (MCI)
% Mixed Effects Inference 
% FORMAT [MCI] = spm_mci_mfx (MCI)
%
% MCI               Data structure containing fields:
%
% .M{n}             Model for nth of N replications (e.g. subjects)
% .U{n}             Inputs for nth replication
% .Y{n}             Data for nth replication
% .S                Second level model describing population mean, m, and 
%                   precision, Lambda. The parameters in S.prior
%                   define the sufficient statistics of p(Lambda) (.a and .B)
%                   and p(m|Lambda) (.beta and.m)
%
% .inference        'amc' or 'lgv' (default)
% .total_its        Total number of samples per subject
% .rinit            Proportion of samples to collect prior to use of
%                   Empirical (group) prior
% .verbose          Show progress of optimisation 
% .update_obs_noise Update observation noise ? [yes/no] (1/0), default=1
%
% The output fields are: 
%
% POSTERIOR SAMPLES:
% .sm               [Nw x Nsamples] group random effect means, m
% .sw               [Nw x N x Nsamples] subject random effects, w
% .Ce               [Ny x Ny x Nsamples] Obs noise covariance samples
% .postind          Indices for posterior (ie. excluding burn-in)
%
% POSTERIOR MEANS:
% .sm_mean          [Nw x 1] posterior mean over m
% .sw_mean          [Nw x N] posterior mean over w
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
% $Id: spm_mci_mfx.m 6697 2016-01-27 14:57:28Z spm $

try, verbose=MCI.verbose; catch, verbose=0; end
try, inference=MCI.inference; catch, inference='lgv'; end

% Update observation noise ?
try, update_obs_noise=MCI.update_obs_noise; catch, update_obs_noise=1; end

try, total_its=MCI.total_its; catch, total_its=1024; end
try, rinit=MCI.rinit; catch, rinit=0.25; end

rfx_its=2;
mfx_its=ceil((1-rinit)*total_its);
rfx1_its=total_its-mfx_its;

M=MCI.M;U=MCI.U;Y=MCI.Y;
Np=length(spm_vec(M{1}.pE));

% Observation noise
Ny=size(Y{1}.y,2);
noise.c0=Ny;
noise.D0=eye(Ny);
   
% Prior over random effects (second level model)
try
    S=MCI.S;
catch
    S.prior.P=Np;
    S.prior.a=Np/2;
    S.prior.B=eye(Np);
    S.prior.beta=1;
    S.prior.m=zeros(Np,1);
    S.N=length(M);
end
% Prior mean and cov of random effects from second level model
m=S.prior.m;
C=S.prior.B/S.prior.a;
rfx_prior.pE=m;
rfx_prior.pC=C;
MCI.S=S;

% Initialisation of each model
subj_mcmc.verbose=verbose;
subj_mcmc.maxits=rfx1_its;
for n=1:S.N,
    subj_mcmc.init=M{n}.pE;
    if rfx1_its > 0
        switch inference,
            case 'lgv',
                Psamp = spm_mci_random (subj_mcmc,rfx_prior,[],M{n},U{n},Y{n});
                
            case 'amc',
                M{n}.pE=rfx_prior.pE;
                M{n}.pC=rfx_prior.pC;
                Nscale=0;
                Nsamp=ceil(0.5*rfx1_its);
                Ntune=Nsamp;
                amc = spm_mci_popdef (Nscale,Ntune,Nsamp);
                amc.verbose=verbose;
                Mc{1}=M{n};Uc{1}=U{n};
                Ptmp=spm_mci_pop (amc,Mc,Uc,Y{n});
                Psamp=Ptmp{1}.theta;
                
            otherwise
                disp('Unknown inference method in spm_mci_mfx.m');
        end
        sw(n,:,:) = Psamp;
        w(:,n) = Psamp(:,end);
    else
        sw(n,:,1) = M{n}.pE;
        w(:,n) = M{n}.pE;
    end
end
sw=permute(sw,[2 1 3]);
sm=squeeze(mean(sw,2));

if isfield(M{1},'Ce')
    gauss_noise=1;
    Ce(:,:,1:rfx1_its)=M{1}.Ce;
else
    gauss_noise=0;
end
subj_mcmc.maxits=rfx_its;


% Main Loop
for it=1:mfx_its,
    
    if verbose
        disp(sprintf('RFX iteration %d',it));
    end
        
    if it>1   
        % Update second level params
        S = spm_nwpost (S,w);
        [m,Lambda,C] = spm_nwrnd (S.post,1);
    end
    
    % Updated prior on random effects
    rfx_prior.pE=m;
    rfx_prior.pC=C;
    
    % 1. Update estimates of subject random effects
    if verbose, disp('Updating random effects'); end
    
    for n=1:S.N,
        % Start sampling where left off
        subj_mcmc.init=w(:,n);
        
        switch inference,
            case 'lgv',
                Psamp = spm_mci_random (subj_mcmc,rfx_prior,[],M{n},U{n},Y{n});
                
            case 'amc',
                M{n}.pE=rfx_prior.pE;
                M{n}.pC=rfx_prior.pC;
                Nscale=0; Ntune=0; Nsamp=rfx_its;
                amc = spm_mci_popdef (Nscale,Ntune,Nsamp);
                amc.verbose=verbose;
                Mc{1}=M{n};Uc{1}=U{n};
                Ptmp=spm_mci_pop (amc,Mc,Uc,Y{n});
                Psamp=Ptmp{1}.theta;
                
            otherwise
                disp('Unknown inference method in spm_mci_mfx.m');
        end
        w(:,n) = Psamp(:,end);
    end
    
    % Update observation noise precision
    if update_obs_noise && gauss_noise
        if verbose, disp('Updating observation noise precision'); end
        [noise,M] = spm_mci_obsnoise (w,[],[],noise,M,U,Y);
    end
    if gauss_noise, Ce(:,:,rfx1_its+it)=M{1}.Ce; end
        
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

% Return posterior means
MCI.sm_mean=mean(sm(:,post_ind),2);
MCI.sw_mean=mean(sw(:,:,post_ind),3);
MCI.post_ind=post_ind;

MCI.M=M;
MCI.S=S;
if gauss_noise
    MCI.noise=noise;
    MCI.Ce=Ce; 
end
