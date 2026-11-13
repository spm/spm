function DCM = dts_meeg_csd_lfp_template(model, data_file, cfg)
% Build a four-source CSD LFP DCM for SPM's LFP_BG dataset.
% Authored by Pranay Yadav in 2026

if nargin < 3 || isempty(cfg)
    cfg = struct();
end
if ~isfield(cfg, 'Nmax')
    if isfield(cfg, 'basic_fit_Nmax')
        cfg.Nmax = cfg.basic_fit_Nmax;
    else
        cfg.Nmax = 4;
    end
end

if exist(data_file, 'file') ~= 2
    error('Missing SPM CSD LFP data fixture: %s', data_file);
end

DCM = [];
DCM.xY.Dfile = data_file;
DCM.xY.modality = 'LFP';
DCM.xY.Ic = [2 12 13 18];

DCM.options.analysis = 'CSD';
DCM.options.model = model;
DCM.options.spatial = 'LFP';
DCM.options.trials = 1;
DCM.options.Tdcm = [1 10000];
DCM.options.Fdcm = [2 64];
DCM.options.Rft = 8;
DCM.options.onset = 60;
DCM.options.Nmodes = 4;
DCM.options.h = 1;
DCM.options.han = 0;
DCM.options.D = 2;
DCM.options.lock = 0;
DCM.options.location = 0;
DCM.options.symmetry = 0;
DCM.options.dur = 16;
DCM.options.DATA = 1;
DCM.options.Nmax = cfg.Nmax;

DCM.Sname = {'GP1'; 'GP2'; 'STN1'; 'STN2'};
DCM.Lpos = sparse(3, 0);

DCM.A = {sparse(4, 4), sparse(4, 4), sparse(4, 4)};
DCM.A{1}(3:4, 1:2) = 1; % GP to STN forward connections
DCM.A{2}(1:2, 3:4) = 1; % STN to GP backward connections
DCM.B = {};
DCM.C = sparse(4, 0);

DCM.xU.X = sparse(1, 0);
DCM.xU.name = {};

DCM.M.Nmax = cfg.Nmax;
DCM.M.nograph = 1;
DCM.M.noprint = 1;
end
