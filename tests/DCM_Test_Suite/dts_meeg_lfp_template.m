function DCM = dts_meeg_lfp_template(model, dfile, cfg)
% Base one-node LFP DCM for ERP-family model tests.
% Authored by Pranay Yadav in 2026

if nargin < 2
    dfile = '';
end
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
if ~isfield(cfg, 'integrator'), cfg.integrator = 'spm_int_L'; end

DCM = [];
DCM.xY.Dfile    = dfile;
DCM.xY.modality = 'LFP';

DCM.options.trials   = 1;
DCM.options.Tdcm     = [0 500];
DCM.options.D        = 1;
DCM.options.h        = 0;
DCM.options.han      = 1;
DCM.options.Nmodes   = 1;
DCM.options.analysis = 'ERP';
DCM.options.model    = model;
DCM.options.spatial  = 'LFP';
DCM.options.onset    = 60;
DCM.options.dur      = 16;
DCM.options.CVA      = 0;
DCM.options.Nmax     = cfg.Nmax;
DCM.options.DATA     = 1;
DCM.options.lock     = 0;
DCM.options.multiC   = 0;
DCM.options.symmetry = 0;
DCM.options.location = 0;

DCM.M.Nmax       = cfg.Nmax;
DCM.M.integrator = cfg.integrator;
DCM.M.hE         = 6;
DCM.M.hC         = 1/128;
DCM.M.nograph    = 1;
DCM.M.noprint    = 1;

DCM.Sname = {'A1'};
DCM.Lpos  = sparse(3, 0);
DCM.A     = {1, 1, 0};
DCM.B     = {1};
DCM.C     = 1;

DCM.xU.X    = sparse(1, 0);
DCM.xU.name = {};
end
