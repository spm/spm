% analyse some ERP data (mismatch negativity ERP SPM file from SPM-webpages)
% This is an example batch script to analyse two evoked responses with an
% assumed 5 sources.
% To try this out on your data (the date of this example don't exist in your SPM8 distribution), 
% you have to change 'Pbase' to your own analysis-directory, and choose a name ('DCM.xY.Dfile') 
% of an existing SPM for M/EEG-file with at least two evoked responses. 

% Please replace filenames etc. by your own.
%--------------------------------------------------------------------------
spm('defaults','EEG');

% Data and analysis directories
%--------------------------------------------------------------------------

Pbase     = '.';        % directory with your data, 

Pdata     = fullfile(Pbase, '.'); % data directory in Pbase
Panalysis = fullfile(Pbase, '.'); % analysis directory in Pbase

% Data filename
%--------------------------------------------------------------------------
DCM.xY.Dfile = 'maeMdfspm8_subject1';

% Parameters and options used for setting up model
%--------------------------------------------------------------------------
DCM.options.analysis = 'ERP'; % analyze evoked responses
DCM.options.model    = 'ERP'; % ERP model
DCM.options.spatial  = 'IMG'; % spatial model
DCM.options.trials   = [1 2]; % index of ERPs within ERP/ERF file
DCM.options.Tdcm(1)  = 0;     % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = 200;   % end of peri-stimulus time to be modelled
DCM.options.Nmodes   = 8;     % nr of modes for data selection
DCM.options.h        = 1;     % nr of DCT components
DCM.options.onset    = 60;    % selection of onset (prior mean)
DCM.options.D        = 1;     % downsampling

%--------------------------------------------------------------------------
% Data and spatial model
%--------------------------------------------------------------------------
DCM  = spm_dcm_erp_data(DCM);

%--------------------------------------------------------------------------
% Location priors for dipoles
%--------------------------------------------------------------------------
DCM.Lpos  = [[-42; -22; 7] [46; -14; 8] [-61; -32; 8] [59; -25; 8] [46; 20; 8]];
DCM.Sname = {'left AI', 'right A1', 'left STG', 'right STG', 'right IFG'};
Nareas    = size(DCM.Lpos,2);

%--------------------------------------------------------------------------
% Spatial model
%--------------------------------------------------------------------------
DCM = spm_dcm_erp_dipfit(DCM);

%--------------------------------------------------------------------------
% Specify connectivity model
%--------------------------------------------------------------------------
cd(Panalysis)

DCM.A{1} = zeros(Nareas,Nareas);
DCM.A{1} = zeros(Nareas, Nareas);
DCM.A{1}(3,1) = 1;
DCM.A{1}(4,2) = 1;
DCM.A{1}(5,4) = 1;

DCM.A{2} = zeros(Nareas,Nareas);
DCM.A{2}(1,3) = 1;
DCM.A{2}(2,4) = 1;
DCM.A{2}(4,5) = 1;

DCM.A{3} = zeros(Nareas,Nareas);
DCM.A{3}(4,3) = 1;
DCM.A{3}(3,4) = 1;

DCM.B{1} = DCM.A{1} + DCM.A{2};
DCM.B{1}(1,1) = 1;
DCM.B{1}(2,2) = 1;

DCM.C = [1; 1; 0; 0; 0];

%--------------------------------------------------------------------------
% Between trial effects
%--------------------------------------------------------------------------
DCM.xU.X = [0; 1];
DCM.xU.name = {'rare'};

%--------------------------------------------------------------------------
% Invert
%--------------------------------------------------------------------------
DCM.name = 'DCMexample';

DCM      = spm_dcm_erp(DCM);
