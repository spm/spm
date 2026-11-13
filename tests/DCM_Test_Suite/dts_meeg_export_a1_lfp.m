function outputs = dts_meeg_export_a1_lfp(dcm_file, output_file, opts)
% Make a single bilateral A1 LFP-like M/EEG dataset from a fitted ERP DCM.
% Authored by Pranay Yadav in 2026
%
% Inputs:
%   dcm_file     - path to a fitted ERP DCM .mat
%   output_file  - path to write the LFP M/EEG .mat
%   opts (struct, optional fields and defaults):
%       a1_source_indices    : [1 2]    indices of lA1/rA1 in DCM.Sname
%       Nm                   : 1        local source modes per ROI
%       alpha                : 0.005    ridge regularization
%       anchor_first_sample  : true     subtract sample-1 from each LFP
%       overwrite_outputs    : true     remove existing output before writing
%
% Output struct:
%       lfp_file         - path to written LFP M/EEG .mat
%       diagnostics_file - path to diagnostics .mat
%
% Pipeline:
%   1. Load fitted ERP DCM,
%   2. Rebuild IMG regional lead field for lA1 and rA1 only,
%   3. Estimate lA1/rA1 source waveforms from DCM-prepared sensor data using
%      a ridge-regularized minimum-norm inverse,
%   4. Sign-align the two A1 waveforms,
%   5. Average them into one waveform per condition,
%   6. Write the waveform as an SPM M/EEG object with channel type LFP.

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'a1_source_indices'),   opts.a1_source_indices = [1 2]; end
if ~isfield(opts, 'Nm'),                  opts.Nm = 1; end
if ~isfield(opts, 'alpha'),               opts.alpha = 0.005; end
if ~isfield(opts, 'anchor_first_sample'), opts.anchor_first_sample = true; end
if ~isfield(opts, 'overwrite_outputs'),   opts.overwrite_outputs = true; end

a1_source_indices   = opts.a1_source_indices;
Nm                  = opts.Nm;
alpha               = opts.alpha;
anchor_first_sample = opts.anchor_first_sample;
overwrite_outputs   = opts.overwrite_outputs;

%% SPM setup
spm('defaults', 'eeg');
spm_get_defaults('cmdline', true);
spm('CmdLine', true);

%% Load the fitted DCM
S = load(dcm_file, 'DCM');
DCM = S.DCM;

out_dir = fileparts(output_file);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

if overwrite_outputs
    if exist(output_file, 'file')
        Dold = spm_eeg_load(output_file);
        delete(Dold);
    end
    diagnostics_file = spm_file(output_file, 'suffix', '_diagnostics');
    if exist(diagnostics_file, 'file')
        delete(diagnostics_file);
    end
end

fprintf('Input DCM:\n%s\n', dcm_file);
fprintf('Output LFP M/EEG file:\n%s\n', output_file);

%% Rebuild IMG regional lead fields with one mode per A1 ROI
DCM.options.spatial = 'IMG';
if ~isfield(DCM, 'M') || ~isstruct(DCM.M)
    DCM.M = struct();
end
if ~isfield(DCM.M, 'dipfit') || ~isstruct(DCM.M.dipfit)
    DCM.M.dipfit = struct();
end
DCM.M.dipfit.Nm = Nm;
DCM = spm_dcm_erp_dipfit(DCM, 1);

roi_leadfields = DCM.M.dipfit.G(a1_source_indices);
CL = full(cat(2, roi_leadfields{:}));                  % channels x 2

% Global RMS lead-field normalization
Lscale = sqrt(trace(CL * CL') / size(CL, 1));
SL = CL / Lscale;

lambda = alpha * trace(SL' * SL) / size(SL, 2);
W = SL / (SL' * SL + lambda * eye(size(SL, 2)));  % channels x 2

diagnostics = [];
diagnostics.dcm_file = dcm_file;
diagnostics.a1_source_indices = a1_source_indices;
diagnostics.a1_source_names = DCM.Sname(a1_source_indices);
diagnostics.Nm = Nm;
diagnostics.alpha = alpha;
diagnostics.lambda = lambda;
diagnostics.Lscale = Lscale;
diagnostics.leadfield = CL;
diagnostics.scaled_leadfield = SL;
diagnostics.inverse_operator = W;
diagnostics.leadfield_correlation = corrcoef(SL);
diagnostics.singular_values = svd(SL);
diagnostics.condition_number = cond(SL);
diagnostics.resolution_matrix = SL' * W;
diagnostics.noise_gain = sqrt(diag(W' * W));

%% Invert each DCM-prepared condition to lA1/rA1 and average
source_a1 = cell(size(DCM.xY.y));
lfp = cell(size(DCM.xY.y));
lfp_raw_average = cell(size(DCM.xY.y));
time_seconds = DCM.xY.pst(:)' / 1000;

for c = 1:numel(DCM.xY.y)
    Y = DCM.xY.y{c};                 % time x sensors, DCM-prepared data
    source_a1{c} = Y * W;            % time x 2, columns are lA1/rA1

    if source_a1{c}(:,1)' * source_a1{c}(:,2) < 0
        source_a1{c}(:,2) = -source_a1{c}(:,2);
    end

    lfp_raw_average{c} = mean(source_a1{c}, 2);
    lfp{c} = lfp_raw_average{c};

    if anchor_first_sample
        lfp{c} = lfp{c} - lfp{c}(1);
    end
end

%% Convert the one-channel condition waveforms to an SPM M/EEG LFP object
ftdata = [];
ftdata.fsample = 1 / DCM.xY.dt;
ftdata.label = {'A1_LFP'};
ftdata.label = ftdata.label(:);

for c = 1:numel(lfp)
    ftdata.trial{c} = lfp{c}(:)';     % FieldTrip expects channels x time.
    ftdata.time{c} = time_seconds;
end

Dlfp = spm_eeg_ft2spm(ftdata, output_file);
Dlfp = type(Dlfp, 'evoked');
Dlfp = chantype(Dlfp, ':', 'LFP');
Dlfp = units(Dlfp, ':', 'a.u.');

if isfield(DCM.xY, 'code') && numel(DCM.xY.code) == numel(lfp)
    for c = 1:numel(lfp)
        Dlfp = conditions(Dlfp, c, char(DCM.xY.code{c}));
    end
else
    for c = 1:numel(lfp)
        Dlfp = conditions(Dlfp, c, sprintf('condition_%d', c));
    end
end

Dlfp = history(Dlfp, 'dts_meeg_export_a1_lfp', struct( ...
    'dcm_file', dcm_file, ...
    'a1_source_indices', a1_source_indices, ...
    'Nm', Nm, ...
    'alpha', alpha, ...
    'lambda', lambda, ...
    'anchor_first_sample', anchor_first_sample));
save(Dlfp);

diagnostics_file = spm_file(output_file, 'suffix', '_diagnostics');
save(diagnostics_file, ...
    'diagnostics', 'source_a1', 'lfp_raw_average', 'lfp', 'time_seconds');

fprintf('\nFinished A1 LFP export.\n');
fprintf('M/EEG file:    %s\n', fullfile(Dlfp.path, Dlfp.fname));
fprintf('Diagnostics:   %s\n', diagnostics_file);

outputs = struct();
outputs.lfp_file = fullfile(Dlfp.path, Dlfp.fname);
outputs.diagnostics_file = diagnostics_file;
end
