% Example code for processing the data and training the model
% Step 1. Run the MB fitting to generate warps (y_*.nii)
%         and model parameters.
% Step 2. Use the model parameters to generate tissue maps
%         for each training image (c01*.nii, c02*.nii and c03*.nii).
% Step 3. Convert the tissue maps and labels into a form suitable
%         for training the model (pcat_*.nii and pcat_*.mat).
% Step 4. Train the model.
%================================================================================


% Add required SPM12 toolboxes to the search path
toolboxes = {'mb','Shoot','Longitudinal'};
tbx = fullfile(spm('dir'),'toolbox');
for i=1:numel(toolboxes)
    path(fullfile(tbx,toolboxes{i}),path);
end

% Directory containing "images" and "labels" subdirectories
datadir  = '.';
labeldir = fullfile(datadir,'labels');    % Label maps
imagedir = fullfile(datadir,'images');    % Images
tisdir   = fullfile(datadir,'processed'); % MB output

%================================================================================
% STEP 1
% Set up a job to run the MB (MultiBrain) registration.
% Requires files generated previously by the MB toolbox:
%   mu_X.nii      - Template
%                   From: https://figshare.com/s/32dc759911b3ec4a1f11
%   prior_X_3.mat - priors on intensity distribution
%                   From https://figshare.com/s/fda9d680fee1a6fe1314
clear matlabbatch
clear chan
chan.images      = cellstr(spm_select('FPlist',imagedir,'^.*\.nii')); % Images
chan.inu.inu_reg = 10000;               % Intensity non-uniformity (INU) regularisation
chan.inu.inu_co  = 40;                  % Cut-off for INU
chan.modality    = 1;                   % Indicates MRI

clear gmm
gmm.chan           = chan;              % Attach
gmm.labels.false   = [];                % Not informed by labels
gmm.pr.file        = {'prior_X_3.mat'}; % Uses intensity priors
gmm.pr.hyperpriors = [];                % Do not update priors
gmm.tol_gmm        = 0.0005;            % When to stop clustering (GMM)
gmm.nit_gmm_miss   = 32;                % Maximum iterations for missing data
gmm.nit_gmm        = 8;                 % Maximum clustering (GMM) iterations
gmm.nit_appear     = 4;                 % Number of times to iterate the GMM and INU

clear run
run.mu.exist   = {'mu_X.nii'};          % Tissue priors
run.aff        = 'SE(3)';               % Include rigid alignment
run.v_settings = [0.0001 0.5 0.5 0 1];  % Registration regularisation
run.onam       = 'traindat';            % For naming the various outputs
run.odir       = {tisdir};              % Where to write the various outputs
run.cat        = {{}};                  % Not working with pre-segmented images
run.gmm        = gmm;                   % Attach
run.accel      = 0.8;                   % 0 is slow and stable. 1 is fast and unstable
run.min_dim    = 8;                     % Minimum image size for multi-scale fitting
run.tol        = 0.001;                 % When to stop
run.sampdens   = 2;                     % Speed accuracy tradeoff
run.save       = true;                  % Save everything to disk
run.nworker    = 0;                     % For parallelisation
matlabbatch{1}.spm.tools.mb{1}.run = run; % Attach
%================================================================================

%================================================================================
% STEP 2.
% Settings to generate c01*.nii, c02*.nii and c03*.nii tissue maps
clear out
out.result  = {fullfile(run.odir{1},['mb_fit_' run.onam '.mat'])};
out.i       = false;                    % Don't generate images (with missing data filled in)
out.mi      = false;                    % Don't generate INU corrected images
mb.out.wi   = false;                    % Don't generate warped and INU corrected images
out.wmi     = false;                    % Don't generate modulated, warped and INU corrected images
out.inu     = false;                    % Don't generate INU maps
out.c       = [1 2 3];                  % Save tissues 1, 2 and 3
out.wc      = [];                       % Don't save warped tissues
out.mwc     = [];                       % Don't save modulated and warped tissues
out.sm      = [];                       % Don't save scalar momenta maps
out.mrf     = 0;                        % Don't use MRF corrections
out.fwhm    = [0 0 0];                  % Don't smooth
out.bb      = [NaN NaN NaN
               NaN NaN NaN];            % Use default bounding box for any warped images
out.vox     = NaN;                      % Use default voxel sizes for any warped images
out.proc_zn = {};
out.odir    = {run.odir{1}};            % Output directory
matlabbatch{2}.spm.tools.mb.out = out;  % Attach

% Run the job
spm_jobman('run',matlabbatch)
%================================================================================

%================================================================================
% STEP 3.
% Generate the pcat_* from the c01*.nii, c02*.nii, c03*.nii, labels and y_*.nii.

% Get a suitable bounding box
mb      = out.result{1};
Pmu     = mb.sett.mu.exist.mu;
[dw,Mw] = fil_subvol(Pmu,[28 29 43; 170 209 194]);

% Get warps
Pdef    = spm_select('FPList',tisdir,['^y_1_.*\.nii$']);
Niiy    = nifti(Pdef);
Niiy    = Niiy(:);

% Generate pushed tissue maps
Ntis    = [];
for c=1:3
    Ptis   = spm_select('FPList',tisdir,['^c0' num2str(c) '.*\.nii$']);
    Ntc    = nifti(Ptis);
    Ntis   = [Ntis Ntc(:)];
end
fil_push_train_data(dw, Mw, Niiy, Ntis);

% Generate pushed label maps
Plabel   = spm_select('FPList',labeldir,'^1.*_glm\.nii');
Nlabel   = nifti(Plabel);
Nlabel   = Nlabel(:);
fil_push_train_data(dw, Mw, Niiy, Nlabel);
%================================================================================

%================================================================================
% Step 4.
% Train the model using the pcat_*.
% Note that it also generates a bunch of fil_*latent*.mat, which are used during
% training, but not needed at deployment time. 

% Note that this same function was used to train from the MICCAI dataset, but a
% slight change was made to the fil_train.m function so that repeat subjects only
% had half the weighting of subjects that were not repeated.

% Training data - assumed to be in current directory
Plab  = spm_select('FPList','.','^pcat_.*\.mat$');  % Modulated and warped labels
Pseg  = spm_select('FPList','.','^pcat_.*\.nii$');  % Modulated and warped tissues
files = {Pseg,Plab};

% Run the training. Takes a few days.
nam   = 'fil30-nuNaN-v1-d4-K24-r2-sd1.5'; % Training filename
sett  = struct('K',   24, ...     % Number of components per patch
               'nit',  5, ...     % Number of alternating E-M steps per iteration
               'nu0',NaN, ...     % Regularisation strength
               'v0',   1, ...     % Regularisation variance
               'd1',   4, ...     % Patch size
               'nit0', 4, ...     % Number of outer iterations
               'matname',[nam '.mat'], ... % Results file name
               'workers',0, ...   % For possible parallelisation
               'verb', 0, ...     % Verbosity
               'r',  2.0, ...     % Radius for augmentation translations
               'sd',1.5);         % Gaussian weighting for augmentation
fil_train(files,sett);
%================================================================================


