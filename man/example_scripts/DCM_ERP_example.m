% analyse some ERP data
% This is an example batch script to analyse a single ERP.

% Please replace filenames etc. by your own. (The example names are not
% part of the current distribution)
%-------------------------------------------------------
clear all

% paths to data, etc.
%-------------------------------------------------------
Pbase = 'D:\DCM'; % supraordinate directory 

Pdata = fullfile(Pbase, '.'); % data directory in Pbase

Panalysis = fullfile(Pbase, '.'); % analysis directory in Pbase

Dfile = 'cmafwde_nervestim.mat'; % the data (ERP SPM-matfile
DCM.M.dipfit.sensorfile = fullfile(Pbase, 'olivier_09_03_2004.pol'); % a sensor location file

% Parameters and options used for setting up model.
%-------------------------------------------------------
DCM.options.Y1 = 1; % index of 1st ERP within ERP/ERF file
DCM.options.Y2 = 2; % index of 2nd ERP within ERP/ERF file
DCM.options.Tdcm(1) = 5; % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2) = 150; % end of peri-stimulus time to be modelled
DCM.options.Nmodes = 3; % nr of modes for data selection
DCM.options.h = 2; % nr of DCT components
DCM.options.Spatial_type = 1; % spatial model is ECD/EEG (2: ECD/MEG, 3: fixed)
DCM.M.onset = 60; % selection of onset (prior mean)
%----------------------------------------------------------
% data and observation level
%----------------------------------------------------------
%  ERP file
D = spm_eeg_ldata(fullfile(Pdata, Dfile));
Ichannels = setdiff(D.channels.eeg, D.channels.Bad);
DCM.M.Ichannels = Ichannels;
% 1st
DCM.Y.xy{1} = squeeze(D.data(Ichannels, :, DCM.options.Y1))';
%DCM.Y.xy{2} = squeeze(D.data(Ichannels, :, DCM.options.Y2))';
DCM.Y.Time = [-D.events.start:D.events.stop]*1000/D.Radc; % ms

% save original data
DCM.Y.xy_original = DCM.Y.xy;
DCM.Y.Time_original = DCM.Y.Time;

% time points to be modelled
[m, T1] = min(abs(DCM.Y.Time-DCM.options.Tdcm(1)));
[m, T2] = min(abs(DCM.Y.Time-DCM.options.Tdcm(2)));
DCM.Y.Time = DCM.Y.Time(T1:T2);
DCM.Y.xy{1} = DCM.Y.xy{1}(T1:T2, :);
%DCM.Y.xy{2} = DCM.Y.xy{2}(T1:T2, :);
DCM.Y.dt = 1000/D.Radc; % ms

% data selection using SVD
P.projection = 1;
P.Nmodes = DCM.options.Nmodes;
DCM      = spm_dcm_eeg_selectdata(DCM, P);

% number of DCT components to model drift at observation level
DCM.Y.h = DCM.options.h;

%----------------------------------------------------------
% specify connectivity model
%----------------------------------------------------------
cd(Panalysis)

Nareas = 3;

DCM.A{1} = zeros(Nareas, Nareas);
DCM.A{1}(2, 1) = 1;

DCM.A{2} = zeros(Nareas,Nareas);
DCM.A{2}(1, 2) = 1;

DCM.A{3} = zeros(Nareas,Nareas);
DCM.A{3}(2,3) = 1;
DCM.A{3}(3,2) = 1;

DCM.B{1} = zeros(Nareas,Nareas);

% DCM.B{1} = DCM.A{1} + DCM.A{2} + DCM.A{3};

DCM.C = [1; 0; 0];

% input 
% DCM.U.X = [0; 1];
DCM.U.X = [1];
DCM.U.name = {'left'};

DCM.Sname = {'left SI', 'left SII', 'right SII'};

% indices of measured channels
DCM.Ichannels = Ichannels;

%----------------------------------------------------------
% specify spatial model
%----------------------------------------------------------

% enter zero leadfield for now, we're using parameterised leadfield 
DCM.L = zeros(DCM.options.Nmodes, Nareas);


% prepare spatial forward model (ECD/EEG) for DCM
DCM = spm_dcm_erp_prepareSpatial(DCM);

dipfit = DCM.M.dipfit;
% transformation matrix from dip to MNI
Mt = [[0 -1 0 0]; [1 0 0 -20]; [0 0 1 -10]; [0 0 0 1]];

% dipole locations for SI/SII areas (MNI coordinates)
dipfit.L.pos = [[-45; -25; 35] [-50; -20; 5] [50; -20; 5]];
dipfit.L.pos = inv(Mt)*[dipfit.L.pos; ones(1, Nareas)];
dipfit.L.pos = dipfit.L.pos(1:3,:);

% priors on dipole orientation as well (zero priors == no orientation
% specified)
dipfit.L.mom = zeros(3, Nareas);

dipfit.L.K = ones(Nareas, 1);

% priors variances
dipfit.L.Vpos = 8*ones(3,Nareas); % tight
dipfit.L.Vmom = 8*ones(3,Nareas); % wide
dipfit.L.VK = ones(size(dipfit.L.K)); % doesn't matter much

DCM.M.dipfit = dipfit;

DCM.M.Lpos = NaN*dipfit.L.pos; % NaN forces spm_gx_erp to compute initial lead field
DCM.M.Lmom = NaN*dipfit.L.mom; % NaN forces spm_gx_erp to compute initial lead field

% necessary for making GUI work
DCM.options.data_ok = 1;
DCM.options.spatial_ok = 1;


DCM.name = '3areas_3comp_batch';

DCM = spm_dcm_erp(DCM);



