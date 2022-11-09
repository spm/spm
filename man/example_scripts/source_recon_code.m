% Example script for running source reconstruction on a single subject from
% Wakeman & Henson dataset as used in SPM MEEG demos.
%__________________________________________________________________________

% Stephanie Mellor
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging


clear

%addpath('**path to spm12***');

%%
spm('defaults', 'eeg');


%% 0. Load Data
%cd('Preprocessing_demo\Sub15');

D = spm_eeg_load(fullfile(pwd, 'maMefdspmeeg_run_01_sss.mat'));

% Remove previous reconstructions
if isfield(D, 'inv')
    D = rmfield(D, 'inv');
end

%% 1.(a) Create cortical, skull and scalp meshes from MRI

sMRI = fullfile(pwd, 'mprage.nii');

% Set reconstruction number
val = 1;

D.val = val;

% Select normal mesh (as opposed to fine = 3 or coarse = 1) 
Msize = 2;

% Create and save mesh
D.inv{val}.mesh = spm_eeg_inv_mesh(sMRI, Msize);

% Set name of inversion
D.inv{val}.comment = {'inversion_2'};

% Display
spm_eeg_inv_checkmeshes(D);

% Save
save(D);

%% 1.(b) Co-register

% Enter MRI fiducial coordinates
nas = [4.10, 111.60, 0.90];
lpa = [-80.40, 20.90, -12.30];
rpa = [78.60, 9.10, -31.00];

mri_fid = D.inv{val}.mesh.fid;
mri_fid.fid.pnt = [nas; lpa; rpa];
mri_fid.fid.label = {'Nasion', 'LPA', 'RPA'}; 

% Get meeg fiducials
meeg_fid = D.fiducials;

% Register
% Get transformation matrix
S = [];
S.sourcefid = meeg_fid;
S.targetfid = mri_fid;
S.useheadshape = 1;
S.template = 0;
M1 = spm_eeg_inv_datareg(S); % M1 = transformation matrix

% Set D.inv{val}.datareg
D.inv{val}.datareg = struct([]);

% Set D.inv{val}.datareg - EEG
D.inv{val}.datareg(1).sensors  = ft_transform_geometry(M1, D.sensors('EEG'));
D.inv{val}.datareg(1).fid_eeg  = ft_transform_geometry(M1, S.sourcefid);
D.inv{val}.datareg(1).fid_mri  = S.targetfid;
D.inv{val}.datareg(1).toMNI    = D.inv{val}.mesh.Affine;
D.inv{val}.datareg(1).fromMNI  = inv(D.inv{val}.datareg(1).toMNI);
D.inv{val}.datareg(1).modality = 'EEG';
    
% Set D.inv{val}.datareg - MEG
D.inv{val}.datareg(2).sensors  = D.sensors('MEG');
D.inv{val}.datareg(2).fid_eeg  = S.sourcefid;
D.inv{val}.datareg(2).fid_mri  = ft_transform_geometry(inv(M1), S.targetfid);
D.inv{val}.datareg(2).toMNI    = D.inv{val}.mesh.Affine*M1;
D.inv{val}.datareg(2).fromMNI  = inv(D.inv{val}.datareg(1).toMNI);
D.inv{val}.datareg(2).modality = 'MEG';

% Display registration
spm_eeg_inv_checkdatareg(D);

%% 1.(c) Compute forward model

% Choose forward models 
D.inv{val}.forward(1).voltype = 'EEG BEM';      % options: {'EEG BEM', '3-Shell Sphere'}
D.inv{val}.forward(2).voltype = 'Single Shell'; % options: {'Single Sphere', 'MEG Local Spheres', 'Single Shell'}

% Compute forward model
D = spm_eeg_inv_forward(D, val);

% Display
spm_eeg_inv_checkforward(D, val);

%% 2. Invert

D.inv{val}.gainmat = 'SPMgainmatrix_maMefdspmeeg_run_01_sss_1.mat';

% Set inversion parameters
inverse.type = 'GS';    % Select Greedy Search
inverse.trials = D.condlist;    % Select all trials
inverse.woi  = [-200 800];  % Select whole window of interest (in ms)
inverse.Han = 1;        % Choose to use Hanning window
inverse.lpf = 0;        % Set high pass filter low point
inverse.hpf = 256;      % Set low pass filter high point
inverse.pQ = {};        % Don't add additional source priors

D.con = 1;
D.inv{val}.inverse  = inverse;

% Select all modalities
[~, list] = modality(D,1,1);    % List all modalities in D
 D.inv{val}.inverse.modality = list;
 
 % Invert
 D = spm_eeg_invert(D);
 
%% 3.(a) Select a contrast window
 
 % Settings
 D.inv{val}.contrast.woi = [100 250];   % Set Time window of interest
 D.inv{val}.contrast.fboi = [10 20];    % Set Frequency band of interest
 D.inv{val}.contrast.type = 'evoked';   % Choose evoked, induced or trials
 
 % Run
 D = spm_eeg_inv_results(D);
 
%% 3.(b) View/save an image

D.inv{val}.contrast.format = 'image';
D = spm_eeg_inv_Mesh2Voxels(D);

save(D);

% Get names of output images
wEEG = D.inv{val}.contrast.fname{D.con};

% Plot
spm_check_registration(fullfile(spm('dir'),'canonical','single_subj_T1.nii')); % Plot template MRI
spm_orthviews('addcolouredimage',1,[wEEG, ',1'],[1 0 0]);   % Overlay contrast image
spm_orthviews('addcolourbar',1,1);                          % Add colorbar
spm_orthviews('Redraw');                                    % Redraw with coloured image                
