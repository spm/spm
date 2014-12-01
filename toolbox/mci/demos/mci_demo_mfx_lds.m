
clear all
close all

disp('Data from multiple subjects');

disp('Estimate mixed effects using Langevin Monte Carlo');
disp('Using LDS model with constrained connectivity');

%lds.model='lds_real';
lds.model='forward';

% Number of states
d=4;

% Observation noise
%lds.sd=0.2;
lds.sd=0.1;

% Prior over initial states
lds.R.pE=linspace(3,1.5,d)';
%lds.R.pC=lds.sd^2*eye(d);
lds.R.pC=0.5^2*eye(d);

% Number of subjects
lds.Nsub=4;

% Number of observations per subject
lds.Nobs=5;

lds.init_par='random';
lds.flow_par='fixed';

% Generate group data
[lds.pinit,lds.pflow,lds.names,M,U,Y] = mci_lds_group_data (lds);

% Assign init/flow as random/fixed effects
assign.init_par='random';
assign.flow_par='fixed';
MCI.assign=assign;

MCI.R=lds.R;

% Initialisation
i0=mci_interp_init(Y,M{1});
a0=[];
% Linear or exponential initialisation ?
% figure;plot(lds.pinit(:),i0(:),'.');
% sum(sum((i0-lds.pinit).^2))
% [i0,a0] = mci_exp_init (Y,M{1});
% figure;plot(lds.pinit(:),i0(:),'.');
% sum(sum((i0-lds.pinit).^2))

if strcmp(assign.init_par,'random')
    MCI.pinit0=i0;
    % Initialise with true values
    %MCI.pinit0=lds.pinit;
else
    MCI.pinit0=mean(i0,2);
    % Initialise with true values
    %MCI.pinit0=mean(lds.pinit,2);
end
if strcmp(assign.flow_par,'random')
    MCI.pflow0=spm_vec(M{1}.pE)*ones(1,lds.Nsub);
    if ~isempty(a0), MCI.pflow0(1:d,:)=a0; end
else
    MCI.pflow0=spm_vec(M{1}.pE);
    if ~isempty(a0), MCI.pflow0(1:d,:)=mean(a0,2); end
end
MCI.verbose=1;
MCI.M=M; MCI.U=U; MCI.Y=Y;

MCI.update_ffx=1;
MCI.update_rfx=1;
MCI.update_obs_noise=1;
MCI.mfx_its=256;

% for n=1:lds.Nsub,
%     lds_plot_fit (MCI,lds,n);
% end

tic;
MCI = spm_mci_mfx (MCI);
MCI = spm_mci_init_flow_xtr (MCI);
toc

rmse=mci_lds_plot_params (MCI,lds);

for n=1:lds.Nsub,
    mci_lds_plot_fit (MCI,lds,n,1);
end

disp('True dynamics:');
[f,Atrue] = mci_lds_fx (lds.pinit(:,1),U{1},lds.pflow,M{1});
disp(Atrue);

disp('Estimated dynamics:');
[f,Aest] = mci_lds_fx (MCI.pinit.Ep,U{1},MCI.pflowK,M{1});
disp(Aest);

if strcmp(assign.flow_par,'fixed')
    disp('Z-scores on flow:');
    pflowK_std=sqrt(diag(MCI.pflowK_cov));
    z=(MCI.pflowK')./pflowK_std
end

% Plot posterior
scale=100;
dist{1}=MCI.pflow;  
dist{1}.color='k';
dist{1}.names=lds.names;
S=size(dist{1}.P,2);
for s=1:S,
    tmp_Pt(:,s)=mci_lds_lat2par (dist{1}.P(:,s),M{1});
end
dist{1}.P=tmp_Pt*scale;
lds_Pt=mci_lds_lat2par(lds.pflow,M{1});
mci_plot_dist_multi (dist,'Flow',lds_Pt*scale);

mci_plot_noiseSD (MCI.Ce,MCI.pflow.ind);
