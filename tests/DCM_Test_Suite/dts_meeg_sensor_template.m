function DCM = dts_meeg_sensor_template(model, dfile, spatial, cfg)
% Base five-source EEG DCM for sensor-level ERP tests.
% Authored by Pranay Yadav in 2026

if nargin < 4 || isempty(cfg)
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
DCM.xY.modality = 'EEG';

DCM.options.trials   = [1 2];
DCM.options.Tdcm     = [0 500];
DCM.options.D        = 1;
DCM.options.h        = 1;
DCM.options.han      = 1;
DCM.options.Nmodes   = 8;
DCM.options.analysis = 'ERP';
DCM.options.model    = model;
DCM.options.spatial  = spatial;
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

DCM.Sname = {'lA1'; 'rA1'; 'lST'; 'rST'; 'rIF'};
DCM.Lpos  = [-42  46 -61  59  46
             -22 -14 -32 -25  20
               7   8   8   8   8];

DCM.A{1} = sparse([0 0 0 0 0
                   0 0 0 0 0
                   1 0 0 0 0
                   0 1 0 0 0
                   0 0 0 1 0]);
DCM.A{2} = sparse([0 0 1 0 0
                   0 0 0 1 0
                   0 0 0 0 0
                   0 0 0 0 1
                   0 0 0 0 0]);
DCM.A{3} = sparse([0 1 0 0 0
                   1 0 1 0 0
                   0 1 0 0 0
                   0 0 0 0 0
                   0 0 0 0 0]);
DCM.B{1} = sparse([1 0 1 0 0
                   0 1 0 1 0
                   1 0 1 0 0
                   0 1 0 1 0
                   0 0 0 0 0]);
DCM.C = sparse([1; 1; 0; 0; 0]);

DCM.xU.X    = [0; 1];
DCM.xU.name = {'deviant'};
end
