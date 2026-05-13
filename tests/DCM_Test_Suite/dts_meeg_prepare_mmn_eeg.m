function outputs = dts_meeg_prepare_mmn_eeg(raw_bdf, raw_pol, final_dir, work_dir, delete_intermediates, overwrite_outputs, make_forward_model)
% Preprocess SPM EEG MMN raw BDF data into an ERP dataset.
% Authored by Pranay Yadav in 2026
%
% Inputs are explicit paths. Outputs is a struct with final_mat, final_dat,
% and, when requested, gainmat_file.

%% Paths and run folder
%--------------------------------------------------------------------------
if nargin < 4
    error('Usage: dts_meeg_prepare_mmn_eeg(raw_bdf, raw_pol, final_dir, work_dir, delete_intermediates, overwrite_outputs)');
end
if nargin < 5 || isempty(delete_intermediates)
    delete_intermediates = true;
end
if nargin < 6 || isempty(overwrite_outputs)
    overwrite_outputs = true;
end
if nargin < 7 || isempty(make_forward_model)
    make_forward_model = false;
end

if overwrite_outputs && exist(final_dir, 'dir')
    rmdir(final_dir, 's');
end
if ~exist(work_dir, 'dir'),  mkdir(work_dir);  end
if ~exist(final_dir, 'dir'), mkdir(final_dir); end

if ~exist(raw_bdf, 'file')
    error('Missing raw BDF file: %s', raw_bdf);
end
if ~exist(raw_pol, 'file')
    error('Missing Polhemus sensor file: %s', raw_pol);
end

fprintf('Raw BDF file:\n%s\n', raw_bdf);
fprintf('Raw sensor file:\n%s\n', raw_pol);
fprintf('Working directory:\n%s\n', work_dir);
fprintf('Final output directory:\n%s\n', final_dir);

%% SPM setup
%--------------------------------------------------------------------------
spm('defaults', 'eeg');
spm_get_defaults('cmdline', true);
spm('CmdLine', true);

%% Raw file hashes
%--------------------------------------------------------------------------
% Deliberately not checked here. The MMN fixture is treated as a standard
% deterministic SPM dataset; if the raw files are wrong, the preprocessing
% steps below will fail visibly.

%% Copy raw inputs into the run folder
%--------------------------------------------------------------------------
copyfile(raw_bdf, fullfile(work_dir, 'subject1.bdf'));
copyfile(raw_pol, fullfile(work_dir, 'sensors.pol'));

raw_bdf = fullfile(work_dir, 'subject1.bdf');
raw_pol = fullfile(work_dir, 'sensors.pol');

cwd = pwd;
cd(work_dir);
cleanup = onCleanup(@() cd(cwd));

%% 1. Convert BioSemi BDF to an SPM continuous M/EEG dataset
% SPM function: spm_eeg_convert
%--------------------------------------------------------------------------
S = [];
S.dataset         = raw_bdf;        % Raw BioSemi BDF file.
S.mode            = 'continuous';   % Read as continuous data.
S.channels        = {'all'};        % Keep EEG plus auxiliary channels for EOG montage.
S.eventpadding    = 0;              % Seconds around events to read additionally.
S.checkboundary   = 1;              % Check for file/data-boundary events.
S.saveorigheader  = 0;              % Do not save the full original header.
S.outfile         = 'spmeeg_subject1';
S.timewin         = [];             % Empty means all available data.
S.conditionlabels = {'Undefined'};  % Continuous data condition label.
S.inputformat     = [];             % Autodetect input format.
D = spm_eeg_convert(S);

%% 2. Load measured EEG sensor positions and fiducials from sensors.pol
% SPM function: spm_eeg_prep
%--------------------------------------------------------------------------
S = [];
S.D        = D.fullfile();
S.task     = 'loadeegsens';
S.source   = 'locfile';
S.sensfile = raw_pol;               % Electrode positions.
D = spm_eeg_prep(S);

%% 3. Build and apply average-reference plus HEOG/VEOG montage
% SPM function: spm_eeg_montage
%--------------------------------------------------------------------------
D = spm_eeg_load(D.fullfile());

montage = [];
montage.labelorg = D.chanlabels(:);
montage.labelnew = [montage.labelorg(1:128), {'HEOG'}, {'VEOG'}];

tra = eye(D.nchannels);
tra(129:end, :) = [];

% Average reference for the 128 EEG channels.
tra = detrend(tra, 'constant');

% HEOG = right peri-orbital minus left peri-orbital in the BioSemi layout.
tra(129, 129) = 0;
tra(129, [131 130]) = [1 -1];

% VEOG = lower eye minus upper eye in the BioSemi layout.
tra(130, 130) = 0;
tra(130, [130 129]) = [1 -1];

montage.tra = tra;
montage.chantypenew = [repmat({'EEG'}, 128, 1); {'EOG'}; {'EOG'}];
montage.chanunitnew = [repmat({'uV'}, 128, 1); {'uV'}; {'uV'}];

S = [];
S.D             = D.fullfile();
S.mode          = 'write';
S.prefix        = 'M';
S.montage       = montage;
S.keepothers    = 0;                % Discard channels not in montage.
S.keepsensors   = 1;                % Transform EEG sensor geometry with montage.
S.updatehistory = 1;
D = spm_eeg_montage(S);

%% 4. High-pass filter at 0.1 Hz
% SPM function: spm_eeg_filter
%--------------------------------------------------------------------------
S = [];
S.D      = D.fullfile();
S.type   = 'butterworth';           % Filter family.
S.band   = 'high';                  % High-pass filter.
S.freq   = 0.1;                     % Cutoff frequency in Hz.
S.dir    = 'twopass';               % Zero-phase forward/backward filtering.
S.order  = 4;                       % Butterworth filter order.
S.prefix = 'f';
D = spm_eeg_filter(S);

%% 5. Low-pass filter at 30 Hz
% SPM function: spm_eeg_filter
%--------------------------------------------------------------------------
S = [];
S.D      = D.fullfile();
S.type   = 'butterworth';           % Filter family.
S.band   = 'low';                   % Low-pass filter.
S.freq   = 30;                      % Cutoff frequency in Hz.
S.dir    = 'twopass';               % Zero-phase forward/backward filtering.
S.order  = 5;                       % Butterworth filter order.
S.prefix = 'f';
D = spm_eeg_filter(S);

%% 6. Epoch standard and oddball trials
% SPM function: spm_eeg_epochs
%--------------------------------------------------------------------------
S = [];
S.D                                  = D.fullfile();
S.timewin                            = [-100 500]; % Peri-stimulus window in ms.
S.trialdef(1).conditionlabel         = 'std';
S.trialdef(1).eventtype              = 'STATUS';
S.trialdef(1).eventvalue             = 1;
S.trialdef(1).trlshift               = 0;
S.trialdef(2).conditionlabel         = 'odd';
S.trialdef(2).eventtype              = 'STATUS';
S.trialdef(2).eventvalue             = 3;
S.trialdef(2).trlshift               = 0;
S.bc                                 = 1;          % Baseline correct after epoching.
S.prefix                             = 'e';
S.eventpadding                       = 0;
D = spm_eeg_epochs(S);

%% 7. Reject artefactual trials using EEG channel threshold
% SPM function: spm_eeg_artefact
%--------------------------------------------------------------------------
S = [];
S.D                          = D.fullfile();
S.mode                       = 'reject';     % Mark bad trials for rejection.
S.badchanthresh              = 0.2;          % Bad-channel proportion threshold.
S.methods.channels           = {'EEG'};      % Apply method to EEG channels.
S.methods.fun                = 'threshchan'; % Per-channel amplitude threshold.
S.methods.settings.threshold = 80;           % Threshold in uV.
S.methods.settings.excwin    = 1000;         % Excision window in ms.
S.append                     = true;         % Append to existing artefact marks.
S.prefix                     = 'a';
D = spm_eeg_artefact(S);

%% 8. Robust average standard and oddball ERPs
% SPM function: spm_eeg_average
%--------------------------------------------------------------------------
S = [];
S.D                    = D.fullfile();
S.robust.ks            = 3;       % Robust averaging tuning constant.
S.robust.bycondition   = false;   % Estimate robust weights across conditions.
S.robust.savew         = false;   % Do not save robust weights.
S.robust.removebad     = true;    % Remove trials marked bad.
S.circularise          = false;   % No circularisation for ERPs.
S.prefix               = 'm';
D = spm_eeg_average(S);

%% 9. Keep only the 128 EEG channels for DCM and source modelling
% SPM function: spm_eeg_montage
%--------------------------------------------------------------------------
D = spm_eeg_load(D.fullfile());
eeg_labels = D.chanlabels(D.indchantype('EEG'))';

montage = [];
montage.labelorg = D.chanlabels(:);
montage.labelnew = eeg_labels(:);
montage.tra      = zeros(numel(eeg_labels), D.nchannels);
for i = 1:numel(eeg_labels)
    montage.tra(i, D.indchannel(eeg_labels{i})) = 1;
end
montage.chantypenew = repmat({'EEG'}, numel(eeg_labels), 1);
montage.chanunitnew = repmat({'uV'},  numel(eeg_labels), 1);

S = [];
S.D             = D.fullfile();
S.mode          = 'write';
S.prefix        = 'eegonly_';
S.montage       = montage;
S.keepothers    = 0;
S.keepsensors   = 1;
S.updatehistory = 1;
D = spm_eeg_montage(S);

%% 10. Write final dataset as EEG_MMN.mat/.dat
% SPM method: copy
%--------------------------------------------------------------------------
D = spm_eeg_load(D.fullfile());
final_file = fullfile(final_dir, 'EEG_MMN.mat');
copy(D, final_file);
D = spm_eeg_load(final_file);

eeg_ind = D.indchantype('EEG');
D = units(D, eeg_ind, 'uV');
save(D);

%% 11. Optionally specify template mesh, EEG BEM forward model, and gain matrix
%--------------------------------------------------------------------------
outputs = [];
outputs.final_mat = final_file;
outputs.final_dat = fullfile(final_dir, 'EEG_MMN.dat');
outputs.gainmat_file = '';

if make_forward_model
    outputs = dts_meeg_add_eeg_bem_forward(final_file);
end

%% 12. Delete intermediate files if requested
%--------------------------------------------------------------------------
if delete_intermediates && exist(work_dir, 'dir')
    rmdir(work_dir, 's');
end

%% 13. Report products
%--------------------------------------------------------------------------
fprintf('\nFinished MMN EEG preprocessing.\n');
fprintf('Final M/EEG file: %s\n', outputs.final_mat);
fprintf('Final data file:  %s\n', outputs.final_dat);
if make_forward_model
    fprintf('Gain matrix:      %s\n', outputs.gainmat_file);
end
end
