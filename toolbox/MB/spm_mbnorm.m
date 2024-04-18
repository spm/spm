function out = spm_mbnorm(job)
% Quick spatial normalisation with MB
% FORMAT spm_mbnorm(P)
% P - an array of filenames of scans (one per subject)
%
% This is intended to show how Multi_brain can be used for
% spatially normalising images.

[mu_file,~,icbm2X] = fil_install;

if isa(job,'struct')
    scans = job.images;
else
    scans = job;
end

% Set search path to find additional functions
if ~isdeployed
    addpath(fullfile(spm('dir'),'toolbox','Shoot'));
    addpath(fullfile(spm('dir'),'toolbox','Longitudinal'));
end

out = struct('def',{cell(length(scans),1)},'warped',{cell(length(scans),1)});
for i=1:length(scans)
    scan          = deblank(scans{i});
    [def,warped]  = do_one(scan, mu_file, icbm2X);
    out.def{i}    = def;
    out.warped{i} = warped;
end


function [def,warped] = do_one(scan, mu_file, icbm2X)
    fprintf('\nWorking on "%s"\n', scan)

    % Run Multi-Brain to align with a pre-existing template
    clear mb
    mb.run.mu.exist = {mu_file};                % Template data
    mb.run.aff = 'SE(3)';                       % Include rigid-body alignment
    mb.run.v_settings = [0.0001 0 0.5 0.2 0.5]; % Deformation regularisation
    mb.run.del_settings = Inf;
    mb.run.onam = 'tmp';                        % For output filenames
    mb.run.odir = {'.'};                        % Output into current directory
    mb.run.cat = {{}};                          % No segmentations to align
    mb.run.gmm.chan.images = {scan};            % Image to align
    mb.run.gmm.chan.inu.inu_reg = 2000;         % INU regularisation
    mb.run.gmm.chan.inu.inu_co = 40; % Cutoff (mm) of INU correction
    mb.run.gmm.chan.modality = 1;    % MRI
    mb.run.gmm.labels.false = [];    % No manual labels
    mb.run.gmm.pr.file = {};         % No priors on the mixture model
    mb.run.gmm.pr.hyperpriors = [];  % No updates of hyperpriors
    mb.run.gmm.tol_gmm = 0.0005;     %  Smaller - more accurate, but slower
    mb.run.gmm.nit_gmm_miss = 1;     % Single channel, so this does not add anything
    mb.run.gmm.nit_gmm = 8;          % Do 8 mixture model iterations per appearance iteration
    mb.run.gmm.nit_appear = 4;       % Do four appearance iterations per warp iteration
    mb.run.accel = 0.8;              % 0.8 seems a good compromise between speed and stability
    mb.run.min_dim = 8;              % Minimum image dimensions to work with over multi-scale
    mb.run.tol = 0.001;              % Smaller - more accurate, but slower
    mb.run.sampdens = 2;             % Another speed/accuracy tradeoff
    mb.run.save = false;             % Don't save the .mat file
    mb.run.nworker = 0;              % No parallelisation

    [odir,onam] = fileparts(mb.run.gmm.chan.images{1}); % Output directory and filename
    [dat,sett]  = spm_mb_init(mb.run);                  % Set up data structure
    dat(1).v    = zeros([dat(1).dm 3],'single');        % Modify structure to work in memory
    dat(1).psi  = zeros([dat(1).dm 3],'single');        % Modify structure to work in memory
    dat(1).onam = onam;
    [dat,sett]  = spm_mb_fit(dat,sett);                 % Run Multi-Brain fitting


    % Convert the contents of dat/sett into what would be saved in a deformation
    % (y_*.nii) file
    M   = sett.mu.Mmu;
    phi = dat(1).psi;
    d   = size(phi);
    phi = reshape(bsxfun(@plus,reshape(phi,[prod(d(1:3)),3])*M(1:3,1:3)',M(1:3,4)'),d);
    mat = spm_dexpm(dat(1).q,sett.B)\sett.mu.Mmu;


    % Invert and compose to get something for spatially normalising
    % to ICBM space, rather than to the average-space image
    [pth,nam,ext] = fileparts(scan);
    %defs.comp{1}.inv.comp{1}.def = {fullfile(pth,['y_1_00001_' nam '_tmp.nii'])};
    defs.comp{1}.inv.comp{1}.supplied.Def = phi;
    defs.comp{1}.inv.comp{1}.supplied.mat = mat;
    defs.comp{1}.inv.space = {mu_file};
    defs.comp{2}.def = {icbm2X};
    defs.out{1}.savedef.ofname = [nam '_icbm.nii']; % Save the composed warp
    defs.out{1}.savedef.savedir.saveusr = {pth};
    defs.out{2}.pull.fnames = {scan};      % Generate a spatially normalised version of the image
    defs.out{2}.pull.savedir.saveusr = {pth};
    defs.out{2}.pull.interp = 1;
    defs.out{2}.pull.mask = 1;
    defs.out{2}.pull.fwhm = [0 0 0];
    defs.out{2}.pull.prefix = '';
    out = spm_deformations(defs);

    def    = out.def{1};
    warped = out.warped{1};


