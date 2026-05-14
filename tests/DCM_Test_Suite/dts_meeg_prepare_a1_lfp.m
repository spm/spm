function lfp_file = dts_meeg_prepare_a1_lfp(eeg_mmn_file, lfp_out_file, work_dir, cfg)
% Derive A1_LFP.mat from a preprocessed EEG_MMN dataset.
% Authored by Pranay Yadav in 2026
%
% Specifies a baseline ERP+IMG DCM on EEG_MMN, prepares its data and
% one-mode IMG lead fields, and then derives the bilateral A1 LFP via
% dts_meeg_export_a1_lfp. No DCM inversion is required for this export.
%
% Inputs:
%   eeg_mmn_file - preprocessed EEG_MMN .mat with forward model + gainmat
%   lfp_out_file - target path for A1_LFP.mat
%   work_dir     - scratch directory for the specified DCM
%   cfg          - cfg-like struct (uses the selected fit Nmax and integrator if present)

if nargin < 4 || isempty(cfg), cfg = struct(); end
if ~isfield(cfg, 'Nmax')
    if isfield(cfg, 'basic_fit_Nmax')
        cfg.Nmax = cfg.basic_fit_Nmax;
    else
        cfg.Nmax = 4;
    end
end
if ~isfield(cfg, 'integrator'), cfg.integrator = 'spm_int_L'; end

if ~exist(eeg_mmn_file, 'file')
    error('dts_meeg_prepare_a1_lfp: missing EEG fixture: %s', eeg_mmn_file);
end
if ~exist(work_dir, 'dir'), mkdir(work_dir); end
out_dir = fileparts(lfp_out_file);
if ~isempty(out_dir) && ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% SPM setup
spm('defaults', 'eeg');
spm_get_defaults('cmdline', true);
spm('CmdLine', true);

%% Stage EEG fixture into work_dir so prep does not mutate the canonical copy
D = spm_eeg_load(eeg_mmn_file);
[~, data_name, data_ext] = fileparts(eeg_mmn_file);
Dwork = copy(D, fullfile(work_dir, [data_name data_ext]));
copy_gainmatrix(D, work_dir);
work_eeg_file = Dwork.fullfile();

%% Build a baseline ERP + IMG DCM and prepare data/lead fields only
DCM = dts_meeg_sensor_template('ERP', work_eeg_file, 'IMG', cfg);
DCM.name = fullfile(work_dir, 'DCM_ERP_IMG_spec_A1_LFP_extraction');
DCM.M.dipfit.Nm = 1;
fprintf('\nPreparing specified ERP+IMG DCM for A1 LFP derivation...\n');
DCM = spm_dcm_erp_data(DCM);
DCM = spm_dcm_erp_dipfit(DCM, 1);
specified_dcm_file = [DCM.name '.mat'];
save(specified_dcm_file, 'DCM', spm_get_defaults('mat.format'));

%% Derive the bilateral A1 LFP
dts_meeg_export_a1_lfp(specified_dcm_file, lfp_out_file, struct('Nm', 1));
lfp_file = lfp_out_file;
end

function copy_gainmatrix(D, out_dir)
try
    gain = D.inv{D.val}.gainmat;
    if isempty(gain), return, end
    [pth, nam, ext] = fileparts(gain);
    if isempty(pth), pth = D.path; end
    src = fullfile(pth, [nam ext]);
    if exist(src, 'file')
        copyfile(src, fullfile(out_dir, [nam ext]));
    end
catch
end
end
