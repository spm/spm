function [isvl] = spm_mci_isvl (mcmc,M,U,Y,VL)
% Compute Log Evidence using Importance Sampling 
% FORMAT [isvl] = spm_mci_isvl (mcmc,M,U,Y,VL)
%
% mcmc          Optimisation parameters  eg.
%
% .maxits       number of samples to use
%
% M             Model structure 
% U             Input structure
% Y             Data 
%
% isvl          
% .logev         log evidence
% .L(s)          log likelihood of sth sample
% .v(s)          importance weight of sth sample
% .logev_est(S)  estimate based on first S samples only
% .logev_boot(b) estimate based on bth bootstrap resample (of size .maxits)
%
% Uses IS with VL posterior as proposal
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_isvl.m 6548 2015-09-11 12:39:47Z will $

M = spm_mci_minit (M);
V  = M.V;
M.vpE=spm_vec(M.pE);

Ep=VL.Ep;
Cp=full(VL.Cp);

VL.post0=-0.5*spm_logdet(Cp);
VL.prior0=-0.5*spm_logdet(M.pC);
VL.iCp=inv(Cp);
VL.ipC=inv(M.pC);

% Generate S samples from prior
S=mcmc.maxits;
w=spm_normrnd(Ep,Cp,S);

parfor s=1:S,
    [L(s),v(s)] = spm_mci_isvl_single (M,U,Y,VL,w(:,s));
end

for s=1:S,
    logev_est(s) = isvl_evidence (L(1:s),v(1:s));
end

% Estimate based on all samples
logev=logev_est(S);

% Normalised importance weights
u = v/sum(v);

% Get logev for bootstrap resamples
Nboot=1000;
for b=1:Nboot,
    ind=ceil(rand(1,S)*S);
    L_resample=L(ind);
    v_resample=v(ind);
    logev_boot(b)=isvl_evidence(L_resample,v_resample);
end

% Outputs
isvl.logev = logev;
isvl.L = L;
isvl.v = v;
isvl.logev_est = logev_est;
isvl.logev_boot = logev_boot;

end

%------------------------------------------------------

function [logev] = isvl_evidence (L,v)

% Normalised importance weights
u = v/sum(v);

% Compute log evidence being careful to avoid overflow
lu = log(u);
lup = lu + L;
lupmax = max (lup);
r = exp(lup - lupmax);
logev = log (sum(r))+lupmax;

end




