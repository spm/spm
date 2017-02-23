% example script of an inversion of a MMC-BGC combined DCM with the generic DCM implementation
% using spectral densities as data feature
% 2 data channels: 1 beamformed MEG source, 1 LFP channel from STN
% 2 conditions: OFF and ON medication
% condition-related modulations allowed in all synaptic connections
% taking posteriors from grand average inversion as mean prior values for G, T, and A
% van Wijk et al. (submitted)

clear all, close all

DCM=[];
DCM.name        = 'example_inversion';
DCM.xY.Dfile    = 'filename';

DCM.Sname       = {'cortex';'basal ganglia'};

DCM.options.model(1).source = 'MMC';    % specify name of the first model (in same order as data)
DCM.options.model(1).B      = 1:14;     % all intrinsic synaptic modulations can be modulated between conditions
DCM.options.model(2).source = 'BGC';    % specify name of the first model (in same order as data)
DCM.options.model(2).J      = 5;        % specify state that contributes to observed data: 1=Striatum, 3=GPe; 5=STN; 7=GPi; 9=Thalamus
DCM.options.model(2).B      = [1:9];    % all intrinsic synaptic modulations can be modulated between conditions

DCM.xY.modality = 'LFP';
DCM.xY.Ic       = 1:2;

DCM.options.Fdcm    = [4 48];
DCM.options.trials  = [1 2];
DCM.options.Tdcm    = [1 3000];
DCM.options.D       = 1;
DCM.options.spatial = 'LFP';

%modulation of individual extrinsic connections - only works if corresponding DCM.Bs are set
DCM.options.model(1).Zf{1,1} = [0 1;0 0];       %single forward modulation
DCM.options.model(1).Zb{1,1} = [0 0;1 0];       %split backward modulation: (in)direct connection
DCM.options.model(1).Zb{2,1} = [0 0;1 0];       %split backward modulation: hyperdirect connection

DCM.A{1}        = zeros(2,2);
DCM.A{1}(1,2)   = 1;                %thalamus -> cortex modeled as forward connection

DCM.A{2}        = zeros(2,2);
DCM.A{2}(2,1)   = 1;                %cortex -> Str / STN modeled as backward connection

DCM.A{3}        = zeros(2,2);

DCM.B{1}        = zeros(2,2);
DCM.B{1}(1,2)   = 1;
DCM.B{1}(2,1)   = 1;

DCM.xU.X        = [0; 1];           %(0 = condition 1, 1 = condition 2)
DCM.xU.name     = {'ONvsOFF'};

DCM.M.nodelay   = 2;                %eliminates all delay self-connections

%%%%%%%%%
%%% OPTIONAL: take priors means from posteriors of inversion grand average
%%% (omit to use defaults prior means)
load('priors_grandaverage_OFF.mat');

DCM.M.MMC_T     = MMC_T;
DCM.M.MMC_G     = MMC_G;
DCM.M.BGC_T     = BGC_T;
DCM.M.BGC_G     = BGC_G;
DCM.M.BGC_E     = BGC_E;
DCM.M.MMC_E     = MMC_E;

%%%%%%%%%%

DCM = spm_dcm_csd(DCM);

%%%%%%%%%
%%% OPTIONAL: use multiple runs to get out of early convergence / local minima
cnt=1;
while cnt<=10;
    
    DCM.M.P=DCM.Ep;
    DCM = spm_dcm_csd(DCM);
    
    cnt=cnt+1;
    disp([fnames{fn},'  rep ',num2str(cnt)]);
end