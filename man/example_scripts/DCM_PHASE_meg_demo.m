
% Fit DCM-phase model to MEG data
%
% With the DCM-GUI one is restricted to models 
% with phase interactions that are sin functions.
% With scripts, this can be relaxed to arbitrary order 
% Fourier series
%
% Additionally with the DCM-GUI one is 
% restricted to only one or two experimental conditions
% This is relaxes in a script

% MEG data file
DCM.xY.Dfile='DIRECTORY\fdas_delay_run1.mat';

NFourier=1;  % Order of Fourier expansion for PIFs

% Conditions 1 and 3; 'control' and 'hard'
DCM.options.trials=[1,3];
DCM.xU.X=[0 1]';

% Just R-MTL, R-VIS, R-IFG
source_pos(1,:)=[27,-18,-27];
source_pos(2,:)=[10,-100,0];
source_pos(3,:)=[39,28,-12];
DCM.Sname{1}='Right-MTL';
DCM.Sname{2}='Right-VIS';
DCM.Sname{3}='Right-IFG';
Nr=length(DCM.Sname);

DCM.Lpos=source_pos';
DCM.options.spatial='ECD';
DCM.xY.modality='MEG';

D=spm_eeg_load(DCM.xY.Dfile);
Ic = D.indchantype(DCM.xY.modality, 'GOOD');
DCM.xY.Ic       = Ic;

% forward model (spatial)
%------------------------------------
DCM = spm_dcm_erp_dipfit(DCM,1);
DCM.M.dipfit.type='ECD';

% Define time points to model 
DCM.options.Tdcm=[5000 6000];

% Project MEG data onto source locations and average
DCM = spm_dcm_phase_data(DCM);

% Define DCM structure
A=[0 1 0; 0 0 0; 0 1 0];
for n=1:NFourier,            
    DCM.As(:,:,n)=A;
    DCM.Ac(:,:,n)=A;
    DCM.Bs{1}(:,:,n)=A;
    DCM.Bc{1}(:,:,n)=A;
end

Ns=size(DCM.xY.y{1},1);

% Define Frequency Range
DCM.options.Fdcm=[4 8];

DCM.xU.u=zeros(Ns,1);

% Fit model
tic;
DCM=spm_dcm_phase(DCM);
toc

% Show results
disp('');
disp('Freqs:');
DCM.Ep.df

disp('Control parameters');
disp('Fitted');
DCM.Ep.As
DCM.Ep.Ac

disp('Memory parameters');
DCM.Ep.As+DCM.Ep.Bs{1}
DCM.Ep.Ac+DCM.Ep.Bc{1}



