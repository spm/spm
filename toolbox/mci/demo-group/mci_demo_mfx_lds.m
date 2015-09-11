
clear all
close all

disp('Data from multiple subjects');

disp('Estimate mixed effects using Langevin Monte Carlo');
disp('Using LDS model with constrained connectivity');

lds.model='forward';

% Number of dynamical states e.g. brain areas
d=4;

% Observation noise
lds.sd=0.1;

% Prior over initial states
lds.R.pE=linspace(3,1.5,d)';
lds.R.pC=0.5^2*eye(d);

% Number of subjects
lds.Nsub=3;

% Number of observations per subject
lds.Nobs=5;

lds.init_par='random';
lds.flow_par='fixed';

% Generate group data
[lds.pinit,lds.pflow,lds.names,M,U,Y] = mci_lds_group_data (lds);

% Assign init/flow as random/fixed effects
assign.init_par='random';
assign.flow_par='fixed';
assign.out_par='known';

MCI.assign=assign;

MCI.fixed.pE=M{1}.pE;
MCI.fixed.pC=M{1}.pC;

% Initialisation
i0=mci_interp_init(Y,M{1});
a0=[];
if strcmp(assign.init_par,'random')
    MCI.pinit0=i0;
else
    MCI.pinit0=mean(i0,2);
end
if strcmp(assign.flow_par,'random')
    MCI.pflow0=spm_vec(M{1}.pE)*ones(1,lds.Nsub);
    if ~isempty(a0), MCI.pflow0(1:d,:)=a0; end
else
    MCI.pflow0=spm_vec(M{1}.pE);
    if ~isempty(a0), MCI.pflow0(1:d,:)=mean(a0,2); end
end
MCI.pout0=[];


MCI.M=M; MCI.U=U; MCI.Y=Y;
MCI.update_obs_noise=1;
MCI.verbose=1;
MCI.total_its=16;
MCI.rinit=0.25;

tic;
MCI = spm_mci_mfx_dynamic (MCI);
toc

rmse=mci_lds_plot_params (MCI,lds);

for n=1:lds.Nsub,
    mci_lds_plot_fit (MCI,lds,n,1);
end

disp('True dynamics:');
[f,Atrue] = mci_lds_fx (lds.pinit(:,1),U{1},lds.pflow,M{1});
disp(Atrue);

disp('Estimated dynamics:');
[f,Aest] = mci_lds_fx (MCI.pinit,U{1},MCI.pflow,M{1});
disp(Aest);

mci_plot_noiseSD (MCI.Ce,MCI.post_ind);
