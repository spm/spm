function [P,L,D] = spm_mci_mh (mcmc,M,U,Y)
% Metropolis Hastings with Gaussian priors and proposals
% FORMAT [P,L,D] = spm_mci_mh (mcmc,M,U,Y)
% 
% mcmc      Optimisation parameters  eg.
%
% .nsamp    number of samples to return 
% .Cprop    proposal density
% .init     initial parameter point
%
% M         Model structure
% U         Inputs
% Y         Data 
%
% P         Posterior samples 
% L         Logjoint history
% D         Diagnostics (D.accept_rate, D.els)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_mh.m 6697 2016-01-27 14:57:28Z spm $

% Compute eigen-parameterisation
M = spm_mci_minit (M);
V  = M.V;

try, init=mcmc.init; catch, init=M.vpE; end

% Initialise proposal covariance
Cprop=V'*mcmc.Cprop*V;
Nr=size(Cprop,1);

% perturbations for proposals
delta_theta=spm_normrnd(zeros(Nr,1),Cprop,mcmc.nsamp);

% Initial param in eigenspace
theta = M.V'*(init-M.vpE);
L = spm_mci_joint (theta,M,U,Y);
                
tic;
accept=1;
for i=2:mcmc.nsamp,
    
    xcand=theta(:,end)+delta_theta(:,i);
    Lnew = spm_mci_joint (xcand,M,U,Y);
    
    % Accept candidate ?
    accept_prob=min(1,exp(Lnew-L(end)));
    if rand(1) < accept_prob,
        theta=[theta xcand];
        L(end+1)=Lnew;
        accept=accept+1;
    else
        theta(:,end+1)=theta(:,end);
        L(end+1)=L(end);
    end

end

% Project parameters from eigenspace back into original space
nj=size(theta,2);
P=M.vpE*ones(1,nj)+V*theta;

D.accept_rate=accept/mcmc.nsamp;
D.els=toc;


