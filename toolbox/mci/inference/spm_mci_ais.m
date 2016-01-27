function [post] = spm_mci_ais (mcmc,M,U,Y,vl)
% Annealed Importance Sampling
% FORMAT [post] = spm_mci_ais (mcmc,M,U,Y,vl)
% 
% mcmc      Optimisation parameters  eg.
%
% .J        number of temperatures
% .anneal   annealing schedule:
%           'sigmoid', 'linear', 'nonlinear', 'log' or 'power'
% .prop     type of proposal: 'lmc' or 'mh' (default)
% .nprop    number of proposals at each temperature
% .maxits   number of independent samples to produce
%
% M         Model structure 
% U         Input structure
% Y         Data 
% vl        Variational Laplace solution
%               .Ep                 Posterior Mean
%               .Cp                 Posterior Covariance
%           If this field is specified then AIS starts sampling
%           from the VL posterior. Otherwise from the model prior.
%
% The function returns data structure 'post' with fields
%
% .P                P(:,j) is jth posterior sample
% .logev            approximation to log evidence
% .logev_se         standard error thereof
% .logev_lower      5th percentile thereof
% .logev_upper      95th percentile thereof
% .logev_resample   resampled log evidences
% .traj             individual trajectories
% .acc              acceptances
% .logw             log of (unnormalised) importance weights
% .q                normalised importance weights
% .E                energy (negative log joint)
% .beta             set of inverse temperatures
%
% R Neal (2001) Annealed Importance Sampling. Statistics and
% Computing, 11, 125-139.
%
% This implementation uses the Matlab Parallel Computing toolbox
% (see use of parfor instead of for below).
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_ais.m 6697 2016-01-27 14:57:28Z spm $
            
try, J=mcmc.J; catch, J=32; end
try, anneal=mcmc.anneal; catch, anneal='power'; end
try, prop=mcmc.prop; catch, prop='mh'; end
try, maxits=mcmc.maxits; catch, maxits=32; end
try, nprop=mcmc.nprop; catch, nprop=5; end
try, verbose=mcmc.verbose; catch, verbose=1; end
try, vl=vl; catch, vl=[]; end;
try, mcmc.rec_traj; catch, mcmc.rec_traj=0; end

nprop=nprop+1;

M = spm_mci_minit (M);
V  = M.V;
Np = size(V,1);
Nsub = size(V,2);

if verbose
    disp(sprintf('Using %s annealing schedule with %d temperatures',anneal,J));
end

% Set annealing schedule
switch anneal
    case 'sigmoid',
        % Good for model switching integration
        x=linspace(-10,10,J);
        beta=1./(1+exp(-x));
    case 'linear',
        beta=linspace(0,1,J);
    case 'geometric',
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
    otherwise
        disp('Unknown type of annealing schedule in spm_mci_ais.m');
        return
end

mcmc.scale=fliplr(10.^-[0:1:nprop-1]);
mcmc.beta=beta;
mcmc.nprop=nprop;

if ~isempty(vl)
    vl.Cp=full(vl.Cp);
    vl.Lambdap=inv(vl.Cp);
    vl.const=-0.5*Np*log(2*pi)-0.5*spm_logdet(vl.Cp);
    % Mean and cov in subspace:
    vl.mr=M.V'*vl.Ep-M.vpE
    vl.Cr=M.V'*vl.Cp*M.V;
end

parfor i=1:maxits,
    if verbose
        disp(sprintf('Acquiring %d out of %d IID samples',i,maxits));
    end
    
    if ~isempty(vl)
        [P(:,i),E(i),logw(i),acc(i,:),traj(i,:,:)] = spm_mci_ais_single_vl (mcmc,M,U,Y,vl);
    else
        [P(:,i),E(i),logw(i),acc(i,:),traj(i,:,:)] = spm_mci_ais_single (mcmc,M,U,Y);
    end
end

% Model evidence
[logev,q] = ais_evidence (logw);

% Now get standard error and lower/upper 90% 
% confidence interval from bootstrapping 
Nboot=1000;
for b=1:Nboot,
    ind=ceil(rand(1,maxits)*maxits);
    logw_resample=logw(ind);
    logev_resample(b)=ais_evidence(logw_resample);
end
logev_se=std(logev_resample);
lind=ceil(0.05*Nboot);
uind=floor(0.95*Nboot);
lsort=sort(logev_resample);
logev_low=lsort(lind);
logev_high=lsort(uind);

% Project parameters from eigenspace of prior back into original space
nj=size(P,2);
P=spm_vec(M.pE)*ones(1,nj)+V*P;

% Get Posterior Mean using importance weights
post.wEp=P*q';

post.P=P;
post.logev=logev;
post.logev_se=logev_se;
post.logev_lower=logev_low;
post.logev_upper=logev_high;
post.logev_resample=logev_resample;
post.acc=acc;
post.traj=traj;
post.q=q;
post.E=E;
post.beta=beta;
post.ind=[1:maxits];
post.logw=logw;

end

%-------------------------------------------------------

function [logev,q] = ais_evidence (logw)

r=max(logw);
u=exp(logw-r);
S=mean(u);
logev=log(S)+r;
q=u/sum(u);

end

