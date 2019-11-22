function out = bf_data
% Prepares the data and initialises the beamforming pipeline
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_data.m 7703 2019-11-22 12:06:29Z guillaume $

% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select a directory where the B.mat file containing the beamforming data will be written.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

D = cfg_files;
D.tag = 'D';
D.name = 'M/EEG dataset';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the M/EEG mat file.'};

val = cfg_entry;
val.tag = 'val';
val.name = 'Inversion index';
val.strtype = 'n';
val.help = {'Index of the cell in D.inv where the forward model is stored.'};
val.val = {1};

gradsource = cfg_menu;
gradsource.tag = 'gradsource';
gradsource.name = 'Where to get MEG sensors';
gradsource.help = {'Taking sensors from D.sensors makes it possible to',...
    'use the same head model and coregistration with multiple datasets.',...
    'This relies on the assumption that the sensors are in head coordinates',...
    'and the fiducals are at the same locations'};
gradsource.labels = {'D.inv', 'D.sensors'};
gradsource.values = {'inv', 'sensors'};
gradsource.val = {'inv'};

space = cfg_menu;
space.tag = 'space';
space.name = 'Coordinate system to work in';
space.help = {'Select the coordinate system for the forward model'};
space.labels = {'MNI-aligned', 'Head', 'Native'};
space.values = {'MNI-aligned', 'Head', 'Native'};
space.val = {'MNI-aligned'};

overwrite = cfg_menu;
overwrite.tag = 'overwrite';
overwrite.name = 'Overwrite BF.mat if exists';
overwrite.help = {'Choose whether to overwrite the existing BF.mat file'};
overwrite.labels = {'Yes', 'No'};
overwrite.values = {1, 0};
overwrite.val = {0};

out = cfg_exbranch;
out.tag = 'data';
out.name = 'Prepare data';
out.val = {dir, D, val, gradsource, space, overwrite};
out.help = {'Prepare the input for beamforming'};
out.prog = @bf_data_run;
out.vout = @bf_data_vout;
out.modality = {'EEG'};
end

function  out = bf_data_run(job)

outdir     = job.dir{1};
val        = job.val;
space      = job.space;
gradsource = job.gradsource;
D      = spm_eeg_load(job.D{1});

if ~isfield(D, 'inv')
    error('Please run head model specification.');
end

if numel(D.inv) < val
    error('Invalid inversion index');
end

% cd(outdir);

%-Ask about overwriting files from previous analyses
%--------------------------------------------------------------------------
if exist(fullfile(outdir, 'BF.mat'),'file') && ~job.overwrite
    str = {'Output directory contains existing BF file:',...
        'Continuing will overwrite existing file!'};
    if spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename);
        fprintf('%-40s: %30s\n\n',...
            'Abort...   (existing BF file)',spm('time'));
        out = []; return
    end
end

BF        = [];

BF.data   = spm_eeg_inv_get_vol_sens(D, val, space, gradsource);

BF.data.D = D;
BF.data.mesh = D.inv{val}.mesh;

try
    delete(fullfile(outdir, 'BF.mat'));
end

bf_save_path(BF,fullfile(outdir, 'BF.mat'));

out.BF{1} = fullfile(outdir, 'BF.mat');

end

function dep = bf_data_vout(job)
% Output is always in field "BF", no matter how job is structured
dep = cfg_dep;
dep.sname = 'BF.mat file';
% reference field "B" from output
dep.src_output = substruct('.','BF');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end
