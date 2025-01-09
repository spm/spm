function DEMO_tvec_erp_mmn(redo)
% DEMO_tvec_erp_mmn: Analyse modulatory dynamics in EEG Mismatch Negativity 
% using Dynamic Causal Modeling of Time-Varying Connectivity.
%
% This script analyzes EEG mismatch negativity (MMN) data using dynamic 
% causal modeling (DCM) in SPM. It includes data preparation, 
% preprocessing, and both initial and advanced DCM analyses to 
% explore neural connectivity and modulatory effects over time.
%
% The script assumes you have access to SPM and the necessary data files. 
% If the data is unavailable, it will download the sample dataset. 
% Customize file paths and parameters as needed.
%
% Key Steps:
% 1. Data Preparation: Prepares directories and checks for necessary files.
% 2. Data Preprocessing: Converts raw EEG data into SPM-compatible format,
% filters, and epochs the data.
% 3. Initial DCM Analysis: Fits a basic DCM model to analyze ERP responses 
% and estimate neural connections.
% 4. Advanced Analysis: Refines the DCM model to assess synaptic plasticity 
% and time-varying connectivity.
% 5. Visualization: Projects connectivity changes over time and compares
% observed vs. modeled ERPs.
%
% Requirements:
% - MATLAB with SPM12 installed.
% - Access to raw EEG data or an internet connection for dataset download.
%
% For further details, refer to the paper:
%   Medrano, J., Friston, K. J., & Zeidman, P. (2024). 
%   Dynamic Causal Models of Time-Varying Connectivity.
%
%__________________________________________________________________________
 
% Johan Medrano
% Copyright (C) 2024-2024 Wellcome Centre for Human Neuroimaging

if nargin == 0, redo = false; end
if isscalar(redo), redo = [redo redo]; end

%% Setup the paths
%==========================================================================

% Please replace filenames etc. by your own.
%--------------------------------------------------------------------------
spm('defaults','EEG');

% Data and analysis directories
%--------------------------------------------------------------------------
Pbase     = '.';                            % directory with your data, 
Pdata     = fullfile(Pbase, 'data');        % data directory in Pbase
Panalysis = fullfile(Pbase, 'analysis');    % analysis directory in Pbase

% Create missing folders
%--------------------------------------------------------------------------
for p = {Pbase Pdata Panalysis}
    if ~exist(p{1}, 'dir')
        mkdir(p{1});
    end
end

% Make paths absolute
%--------------------------------------------------------------------------
Pbase     = getfield(dir(Pbase), 'folder'); 
Pdata     = getfield(dir(Pdata), 'folder'); 
Panalysis = getfield(dir(Panalysis), 'folder'); 


%% Get and prepare data
%==========================================================================

% Move to the data folder
%--------------------------------------------------------------------------
cwd = pwd; 
cd(Pdata); 

% Check for preprocessed data
%--------------------------------------------------------------------------
Dfile = fullfile(Pdata, 'aefdfMspmeeg_subject1.mat'); 

if ~exist(Dfile, "file")

    % If not found, check for source data
    %----------------------------------------------------------------------
    file_path = fullfile(Pdata, 'subject1.bdf'); 

    if ~exist(file_path, 'file')

        % If not found, download dataset
        %------------------------------------------------------------------
        file_path = download_mmn_dataset(Pdata); 
    end
    
    % Preprocess data
    %----------------------------------------------------------------------
    preprocess_mmn_dataset(Pdata, file_path); 
end

% Load preprocessed data 
%--------------------------------------------------------------------------
D = spm_eeg_load(Dfile); 

% Check for spatial model
%--------------------------------------------------------------------------
if ~isfield(D, 'inv')

    % If not there, create it
    %----------------------------------------------------------------------
    D = inverse_mmn_dataset(D.fullfile());
end

cd(cwd);

%% Run a first DCM analysis
%==========================================================================

% Move to the analysis folder
%--------------------------------------------------------------------------
cwd = pwd; 
cd(Panalysis); 

DCM = []; 

% Data filename
%--------------------------------------------------------------------------
DCM.xY.Dfile = Dfile;

% Parameters and options used for setting up model
%--------------------------------------------------------------------------
DCM.options.analysis = 'ERP'; % analyze evoked responses
DCM.options.model    = 'ERP'; % ERP model
DCM.options.spatial  = 'ECD'; % spatial model
DCM.options.trials   = [1 2]; % index of ERPs within ERP/ERF file
DCM.options.Tdcm(1)  = 0;     % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = 400;   % end of peri-stimulus time to be modelled
DCM.options.Nmodes   = 8;     % nr of modes for data selection
DCM.options.h        = 1;     % nr of DCT components
DCM.options.onset    = 60;    % selection of onset (prior mean)
DCM.options.D        = 2;     % downsampling
DCM.options.Nmax     = 128;   % number of E steps

%--------------------------------------------------------------------------
% Data and spatial model
%--------------------------------------------------------------------------
DCM  = spm_dcm_erp_data(DCM);

%--------------------------------------------------------------------------
% Location priors for dipoles
%--------------------------------------------------------------------------
DCM.Lpos  = [[-42; -22; 7] [46; -14; 8] [-61; -32; 8] [59; -25; 8] [46; 20; 8]];
DCM.Sname = {'left A1', 'right A1', 'left STG', 'right STG', 'right IFG'};
Nareas    = size(DCM.Lpos,2);
lA1 = 1; rA1 = 2; lSTG = 3; rSTG = 4; rIFG = 5; 

%--------------------------------------------------------------------------
% Spatial model
%--------------------------------------------------------------------------
DCM = spm_dcm_erp_dipfit(DCM);

%--------------------------------------------------------------------------
% Specify connectivity model
%--------------------------------------------------------------------------
% Forward connections
%--------------------------------------------------------------------------
DCM.A{1} = zeros(Nareas, Nareas);
DCM.A{1}(lSTG, lA1)  = 1;
DCM.A{1}(rSTG, rA1)  = 1;
DCM.A{1}(rIFG, rSTG) = 1;

% Backward connections (opposite the forward connections)
%--------------------------------------------------------------------------
DCM.A{2} = DCM.A{1}';

% Lateral connections
%--------------------------------------------------------------------------
DCM.A{3} = zeros(Nareas,Nareas);
DCM.A{3}(rSTG, lSTG) = 1;
DCM.A{3}(lSTG, rSTG) = 1;

% Connections with modulatory effects
%--------------------------------------------------------------------------
DCM.B{1} = DCM.A{1} + DCM.A{2};
DCM.B{1}(lA1, lA1) = 1;
DCM.B{1}(rA1, rA1) = 1;

% Input
%--------------------------------------------------------------------------
DCM.C = zeros(Nareas, 1);
DCM.C(rA1) = 1; 
DCM.C(lA1) = 1;  

%--------------------------------------------------------------------------
% Between trial effects
%--------------------------------------------------------------------------
DCM.xU.X = [0; 1];
DCM.xU.name = {'rare'};

%--------------------------------------------------------------------------
% Invert
%--------------------------------------------------------------------------
DCM.name = 'DCM_mmn';
if ~exist([DCM.name '.mat'], 'file') || redo(1)
    DCM      = spm_dcm_erp(DCM);
end

cd(cwd); 

%% Find out the repetition number of each tone
%==========================================================================
 
% Load dataset
%--------------------------------------------------------------------------
D = spm_eeg_load(Dfile); 

try
    % Try to load dataset with repetition numbers
    %----------------------------------------------------------------------
    if D.nconditions == 2 
        D = spm_eeg_load(spm_file(D.fullfile(), 'prefix', 'R')); 
    end
catch
    % If fails, compute repetition numbers
    %----------------------------------------------------------------------

    % Find standards
    sti = cellfun(@(x) strcmp(x, 'std'), D.conditions); 
    cn = zeros(1, D.ntrials); 
    for i = 1:D.ntrials
        if sti(i)
            % Maintain a rolling count of the repetition
            if i > 1
                cn(i) = cn(i-1) + 1; 
            else 
                cn(i) = 0; 
            end
        else 
            % Reset on deviant
            if i > 1
                cn(i-1) = -1; 
            end
        end
    end

    % Create a new dataset with repetition numbers as conditions
    %----------------------------------------------------------------------
    Cn = cn; 
    cn = strtrim(cellstr(int2str(cn'))); 
    D = D.conditions(':', cn);
    D = D.fname(spm_file(D.fname, 'prefix', 'C')); 
    save(D); 
end

%% Fit the full DCM with time course of synaptic plasticity 
%==========================================================================

% Move to the analysis folder
%--------------------------------------------------------------------------
cwd = pwd; 
cd(Panalysis); 

% Load the DCM estimated above 
%--------------------------------------------------------------------------
load('DCM_mmn.mat')
DCM0 = DCM; 

% Setup a new DCM
%--------------------------------------------------------------------------
DCM = []; 

% Specify the data  
%--------------------------------------------------------------------------
DCM.xY.Dfile    = fullfile(D); 
DCM.xY.modality = 'EEG';

% Specify options
%--------------------------------------------------------------------------
DCM.options.Tdcm     = DCM0.options.Tdcm; % twin in ms
DCM.options.Nmax     = 128; 
DCM.options.spatial  = 'ECD';     % Equivalent Current Dipole
DCM.options.model    = 'CMC';     % Canonical microcircuit model 
DCM.options.analysis = 'ERP';     % ERP analysis
DCM.options.D        =  2;        % 
DCM.options.Nmodes   = 8;     % nr of modes for data selection
DCM.options.onset    = 60;    % selection of onset (prior mean)

% Specify sources
%--------------------------------------------------------------------------
DCM.Sname = DCM0.Sname;                % Sources names
DCM.Lpos  = DCM0.Lpos;                 % Sources positions, unused

% Specify the conditions to model
%--------------------------------------------------------------------------
[~, I] = sort(D.condlist); 
DCM.options.trials = I;        % Model all conditions

% Get a set of temporal basis functions
%--------------------------------------------------------------------------
bftype = 'Cosine'; 
switch bftype
    case "Fourier"
        bf = []; 
        bf.name = 'Fourier set';
        bf.order= 3; 

    case "Cosine"
        bf = []; 
        bf.name = 'Cosine set';
        bf.order = 6; 

    case "Gamma"
        bf = []; 
        bf.name = 'Gamma functions';
        bf.order = 5; 
end
bf.dt     = 0.5;                         % Setup time interval
bf.length = 0.5*10-bf.dt/2;              % Setup duration
bf        = spm_get_bf(bf);              % Get basis set
bf.bf     = bf.bf(:, 2:end);             % Remove intercept
bf.bf     = bf.bf(1:D.nconditions-1, :); % Remove other conditions
bf.bf     = [zeros(1,size(bf.bf,2)); bf.bf];

% Keep track of basis functions info 
%--------------------------------------------------------------------------
DCM.M.bf = bf; 

% Configure the design matrix
%--------------------------------------------------------------------------
X = bf.bf; 

DCM.xU.X    = X; 
DCM.xU.name = arrayfun(@(i) {sprintf('bf%d', i)}, 1:size(X,2)); 

% Copy the connectivity structure from above
%--------------------------------------------------------------------------
DCM.A = DCM0.A; 
DCM.B = repelem(DCM0.B, size(X, 2)); 
DCM.C = DCM0.C; 

% Prepare the initial parameter set
%--------------------------------------------------------------------------
P   = DCM0.Ep;                                 % Take posterior for mean  
P.B = repelem(DCM0.M.pE.B(1), size(X, 2));     % Fill in modulatory effects 
if strcmp(DCM.options.model, 'CMC') && isfield(DCM0.M.pE, 'N') 
    P.N = repelem(DCM0.M.pE.N(1), size(X, 2)); % For CMC model only
end

% Setup initial parameter estimates
%--------------------------------------------------------------------------
DCM.M.P = P;
DCM.M.Q = DCM0.Eg; 

% Warm up (fill in) DCM fields 
%--------------------------------------------------------------------------
DCM.options.Nmax = 1; 
DCM = spm_dcm_erp(DCM);

% Scale down precision on modulatory effects (defaults are too tight)
%--------------------------------------------------------------------------
for i = 1:size(X, 2)
    DCM.M.pC.B{i} = DCM.M.pC.B{i}*(8); 
end

% Setup number of iterations
%--------------------------------------------------------------------------
DCM.options.Nmax = 128; 
DCM.M.Nmax       = 128; 

%--------------------------------------------------------------------------
% Invert model
%--------------------------------------------------------------------------
DCM.name = 'DCM_mmn_rep';
if ~exist([DCM.name '.mat'], 'file') || redo(2)
    DCM      = spm_dcm_erp(DCM);
end

cd(cwd)

%% Project modulatory effects in the temporal domain
%==========================================================================


% Load the DCM 
%--------------------------------------------------------------------------
load(fullfile(Panalysis, 'DCM_mmn_rep.mat'));

% Extract relevant variables
%--------------------------------------------------------------------------
B           = DCM.B; 
Sname       = DCM.Sname;
M           = DCM.M ;
Qp          = DCM.Ep; 
Cp          = DCM.Cp; 
nconditions = numel(DCM.xY.y);
dt          = DCM.xU.dt; 
pst         = DCM.xY.pst;
X           = DCM.xU.X;
bf          = DCM.M.bf; 

% Upsample the basis set
%--------------------------------------------------------------------------
bfup    = DCM.M.bf; 
bfup.dt = bf.dt /10; 
bfup    = spm_get_bf(bfup); 
Xup     = bfup.bf(:, 2:end); 

bfup.t = (0:bfup.dt:bfup.length)'; 

Xup = Xup - mean(X, 1); 
Xup = Xup ./ (max(Xup, [], 1) - min(Xup, [], 1)); 
Xup = Xup .* (max(X, [], 1) - min(X, [], 1));

bfup.t = bfup.t(1:end-1); 
t = bfup.t / bfup.t(end) * ((nconditions-1)*M.ns*dt + pst(end));

% Project modulatory effects in the time domain 
%--------------------------------------------------------------------------
[Bproj, Cbproj] = spm_dcm_project_B(M, Qp.B, Cp, Xup); 



%% Plot observed and modelled ERPs
%==========================================================================

% Create figure
%--------------------------------------------------------------------------
spm_figure('GetWin', 'Observed and modelled ERPs');clf

ns = length(pst);
dt  = mean(diff(pst)); 

% Plot observed ERPs
%--------------------------------------------------------------------------
subplot(3,1,1)
hold on; 
for i = 1:nconditions
    e = DCM.xY.y{i}; 
    imagesc(pst + (i-1)*ns*dt, 1:size(e, 2), e')
    if i == 1, hold on; end
    xline(64/dt/1e6 + (i-1)*ns*dt, '--r', 'LineWidth', 2);
end 
hold off; axis tight; 
xlabel('Time (s)'); ylabel('Channels'); title('Observed')
legend('', 'Stimulus', 'Location', 'northeastoutside', 'Orientation', 'horizontal')
legend('boxoff')

% generate neuronal response and data
%--------------------------------------------------------------------------
DCM = spm_dcm_erp_dipfit(DCM, 1);              % reload dipfit
L   = feval(DCM.M.G, DCM.Eg,DCM.M);            % get gain matrix
x   = feval(DCM.M.IS,DCM.Ep,DCM.M,DCM.xU);     % prediction (source space)

r = {};
x0  = ones(DCM.M.ns,1)*spm_vec(M.x)';         % expansion point for states
for i = 1:numel(x)
    K{i} = x{i} - x0;                   % centre on expansion point
    r{i} = K{i}*L';                     % signal
end

% Plot modelled ERPs
%--------------------------------------------------------------------------
subplot(3,1,2)
hold on; 
for i = 1:nconditions
    e = r{i}; 
    imagesc(pst + (i-1)*ns*dt, 1:size(e, 2), e')
    if i == 1, hold on; end
    xline(64/dt/1e6 + (i-1)*ns*dt, '--r', 'LineWidth', 2);
end
hold off; axis tight; 
xlabel('Time (s)'); ylabel('Channels'); title('Modelled')
legend('', 'Stimulus', 'Location', 'northeastoutside', 'Orientation', 'horizontal')
legend('boxoff')

%% Plot modulatory effects over time
%==========================================================================

% Create figure
%--------------------------------------------------------------------------
subplot(3,1,3)
legendsitems = {}; 
color = linspecer(12); 

% Plot std. dev. of modulatory trajectory
%--------------------------------------------------------------------------
hold on; 
for i = 1:numel(Bproj)
    low  = Bproj{i} - 1.*sqrt(diag(Cbproj{i,i})); 
    high = Bproj{i} + 1.*sqrt(diag(Cbproj{i,i})); 

    fill([t; flip(t)],[low; flipud(high)], ...
        color(i, :), 'FaceAlpha', 0.15, 'EdgeAlpha', 0);
    legendsitems{i} = ''; 
end

% Plot mean modulatory trajectory
%--------------------------------------------------------------------------
[a,b] =  ind2sub(size(B{1}),find(B{1})); 
for i = 1:numel(Bproj)
    plot(t, Bproj{i}, 'Color', color(i, :), 'LineWidth', 2)
    legendsitems{i+numel(Bproj)} = sprintf('%s to %s', Sname{b(i)}, Sname{a(i)}); 
end

% Setup axis and legend
%--------------------------------------------------------------------------
axis tight; ylim([-2 2]); xlim([0, t(end)]); 
xlabel('Time (s)'); ylabel('Log scaling factor');
legend(legendsitems, 'Location', 'northeastoutside')
legend('boxoff')


end 

%==========================================================================
% Helper functions
%==========================================================================



% download_mmn_dataset
%==========================================================================
function file_path = download_mmn_dataset(download_path)
% Downloads the MMN dataset from SPM website to a target path and return 
% the path to the downloaded file. 
%
%__________________________________________________________________________

if nargin == 0, download_path = 'data'; end

url = "https://www.fil.ion.ucl.ac.uk/spm/download/data/eeg_mmn/subject1.bdf"; 
file_path = fullfile(download_path, 'subject1.bdf');

websave(file_path, url);
end


% preprocess_mmn_dataset
%==========================================================================
function D = preprocess_mmn_dataset(wdir, file_path)
% Preprocesses the MMN dataset. This part is a mere copy from the
% example scripts for that datasets that can be found in the SPM repo
%
%__________________________________________________________________________

spm('defaults', 'eeg');

cwd = pwd;
cd(wdir); 

% Convert to SPM
%--------------------------------------------------------------------------
S = [];
S.dataset = file_path;
S.mode = 'continuous';
S.channels = {'all'};
S.eventpadding = 0;
S.blocksize = 3276800;
S.checkboundary = 1;
S.saveorigheader = 0;
S.outfile = 'spmeeg_subject1';
S.timewin = [];
S.conditionlabels = {'Undefined'};
S.inputformat = [];
D = spm_eeg_convert(S);

% Edit montage
%--------------------------------------------------------------------------
montage.labelorg = D.chanlabels;
montage.labelnew = [montage.labelorg(1:128), 'HEOG', 'VEOG']';

tra = eye(D.nchannels);
tra(129:end, :) = [];

% Create the average reference montage
tra = detrend(tra, 'constant');

% HEOG
tra(129, 129) = 0;
tra(129, [131 130]) = [1 -1];

% VEOG
tra(130, 130) = 0;
tra(130, [130 129]) = [1 -1];

montage.tra = tra;

% Update montage
%--------------------------------------------------------------------------
S = [];
S.D = D.fullfile();
S.mode = 'write';
S.blocksize = 655360;
S.prefix = 'M';
S.montage = montage;
S.keepothers = 0;
S.keepsensors = 1;
S.updatehistory = 1;
D = spm_eeg_montage(S);

% High-pass filter @ 0.1Hz
%--------------------------------------------------------------------------
S = [];
S.D = D.fullfile();
S.type = 'butterworth';
S.band = 'high';
S.freq = 0.1;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';
D = spm_eeg_ffilter(S);

% Downsample @ 200Hz
%--------------------------------------------------------------------------
S = [];
S.D = D.fullfile();
S.fsample_new = 200;
S.prefix = 'd';
D = spm_eeg_downsample(S);

% Low-pass filter @ 30Hz
%--------------------------------------------------------------------------
S = [];
S.D = D.fullfile();
S.type = 'butterworth';
S.band = 'low';
S.freq = 30;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';
D = spm_eeg_ffilter(S);

% Epoch 
%--------------------------------------------------------------------------
S = [];
S.D = D.fullfile();
S.timewin = [-100 400];
S.trialdef(1).conditionlabel = 'std';
S.trialdef(1).eventtype = 'STATUS';
S.trialdef(1).eventvalue = 1;
S.trialdef(1).trlshift = 0;
S.trialdef(2).conditionlabel = 'odd';
S.trialdef(2).eventtype = 'STATUS';
S.trialdef(2).eventvalue = 3;
S.trialdef(2).trlshift = 0;
S.bc = 1;
S.prefix = 'e';
S.eventpadding = 0;
D = spm_eeg_epochs(S);

% Reject bad trials 
%--------------------------------------------------------------------------
S = [];
S.D = D.fullfile();
S.mode = 'reject';
S.badchanthresh = 0.2;
S.methods.channels = {'EEG'};
S.methods.fun = 'threshchan';
S.methods.settings.threshold = 80;
S.methods.settings.excwin = 1000;
S.append = true;
S.prefix = 'a';
D = spm_eeg_artefact(S);

cd(cwd);
end

% inverse_mmn_dataset
%==========================================================================
function D = inverse_mmn_dataset(file_path)
% Specifies the forward model for the dataset and inverts it. 
%
%__________________________________________________________________________


% Compute trial average
%--------------------------------------------------------------------------
S = [];
S.D = file_path;
S.robust.ks = 3;
S.robust.bycondition = false;
S.robust.savew = false;
S.robust.removebad = true;
S.circularise = false;
S.prefix = 'm';
D = spm_eeg_average(S);

% Specify forward model 
%--------------------------------------------------------------------------
matlabbatch{1}.spm.meeg.source.headmodel.D = {D.fullfile()};
matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'spmnas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'spmlpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'spmrpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';

% Invert model
%--------------------------------------------------------------------------
matlabbatch{2}.spm.meeg.source.invert.D(1) = cfg_dep('Head model specification: M/EEG dataset(s) with a forward model', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
matlabbatch{2}.spm.meeg.source.invert.val = 1;
matlabbatch{2}.spm.meeg.source.invert.whatconditions.all = 1;
matlabbatch{2}.spm.meeg.source.invert.isstandard.standard = 1;
matlabbatch{2}.spm.meeg.source.invert.modality = {'All'};

% Run batch
%--------------------------------------------------------------------------
spm_jobman('run', matlabbatch)

% Load epoched and trial-averaged datasets
%--------------------------------------------------------------------------
Dm = spm_eeg_load(D.fullfile());
D  = spm_eeg_load(file_path); 

% Set epoched inverse model from the trial-averaged one
%--------------------------------------------------------------------------
D.inv = Dm.inv; 

% Remove fields from evoked inversion
%--------------------------------------------------------------------------
D.inv{1}.inverse.J = {}; 
D.inv{1}.inverse.trials = {'Undefined'}; 
D.inv{1}.inverse.pst = D.time; 
D.inv{1}.inverse.It  = 1:length(D.time);

% Save
%--------------------------------------------------------------------------
save(D)

end