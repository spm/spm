% analyse some ERP data
% This is an example batch script to analyse a single ERP.

% Please replace filenames etc. by your own.
%-------------------------------------------------------
clear all
spm_defaults

% paths to data, etc.
%-------------------------------------------------------
Pbase     = 'D:\DCM';             % supraordinate directory 
Pdata     = fullfile(Pbase, '.'); % data directory in Pbase
Panalysis = fullfile(Pbase, '.'); % analysis directory in Pbase

% the data (ERP SPM-matfile
DCM.xY.Dfile = 'cmafwde_nervestim.mat';

% a sensor location file
DCM.M.dipfit.sensorfile = fullfile(Pbase, 'olivier_09_03_2004.pol');

% Parameters and options used for setting up model.
%-------------------------------------------------------
DCM.options.type    = 1;      % spatial model is ECD (2: imaging, 3: fixed)
DCM.options.trials  = [1 2];  % index of ERPs within ERP/ERF file
DCM.options.Tdcm(1) = 5;      % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2) = 150;    % end of peri-stimulus time to be modelled
DCM.options.Nmodes  = 3;      % nr of modes for data selection
DCM.options.h       = 0;      % nr of DCT components
DCM.options.onset   = 60;     % selection of onset (prior mean)


%----------------------------------------------------------
% data and spatial model
%----------------------------------------------------------
DCM  = spm_dcm_erp_data(DCM);
DCM  = spm_dcm_erp_dipfit(DCM);

% location priors for dipoles
%----------------------------------------------------------
DCM.M.dipfit.L.pos = [[-45; -25; 35] [-50; -20; 5] [50; -20; 5]];
DCM.Sname = {'left SI', 'left SII', 'right SII'};
Nareas    = size(DCM.M.dipfit.L.pos,2);

%----------------------------------------------------------
% specify connectivity model
%----------------------------------------------------------
cd(Panalysis)

DCM.A{1} = zeros(Nareas,Nareas);
DCM.A{1}(2, 1) = 1;

DCM.A{2} = zeros(Nareas,Nareas);
DCM.A{2}(1, 2) = 1;

DCM.A{3} = zeros(Nareas,Nareas);
DCM.A{3}(2,3) = 1;
DCM.A{3}(3,2) = 1;

DCM.B{1} = zeros(Nareas,Nareas);

DCM.C    = [1; 0; 0];

%----------------------------------------------------------
% between trial effects
%----------------------------------------------------------
DCM.xU.X = [1];
DCM.xU.name = {'left'};

%invert
%----------------------------------------------------------
DCM.name = '3areas_3comp_batch';

DCM      = spm_dcm_erp(DCM);



