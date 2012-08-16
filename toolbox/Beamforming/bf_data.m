function out = bf_data
% Prepares the data and initialises the beamforming pipeline
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_data.m 4847 2012-08-16 17:29:23Z vladimir $

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

space = cfg_menu;
space.tag = 'space';
space.name = 'Coordinate system to work in';
space.help = {'Select the coordinate system for the forward model'};
space.labels = {'MNI-aligned', 'Head', 'Native'};
space.values = {'MNI-aligned', 'Head', 'Native'};
space.val = {'MNI-aligned'};

out = cfg_exbranch;
out.tag = 'bf_data';
out.name = 'Prepare data';
out.val = {dir, D, val, space};
out.help = {'Prepare the input for beamforming'};
out.prog = @bf_data_run;
out.vout = @bf_data_vout;
out.modality = {'EEG'};
end

function  out = bf_data_run(job)

outdir = job.dir{1};
val    = job.val;
space  = job.space;
D      = spm_eeg_load(job.D{1});

if ~isfield(D, 'inv')
    error('Please run head model specification.');
end

if numel(D.inv) < val
    error('Invalid inversion index');
end

cd(outdir);

%-Ask about overwriting files from previous analyses
%--------------------------------------------------------------------------
if exist(fullfile(pwd,'BF.mat'),'file')
    str = {'Current directory contains existing BF file:',...
        'Continuing will overwrite existing file!'};
    if spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename);
        fprintf('%-40s: %30s\n\n',...
            'Abort...   (existing BF file)',spm('time'));
        out = []; return
    end
end

BF        = [];
BF.data.D = D;
BF.data.mesh = D.inv{val}.mesh;

eegind = 0;
megind = 0;
for m = 1:numel(D.inv{val}.forward)
    if strncmp('EEG', D.inv{val}.forward(m).modality, 3)
        eegind = m;
    elseif strncmp('MEG', D.inv{val}.forward(m).modality, 3)
        megind = m;
    end
end

istemplate = D.inv{val}.mesh.template;

if megind > 0
    
    vol      = D.inv{val}.forward(megind).vol;
    datareg  = D.inv{val}.datareg(megind);
    sens     = datareg.sensors;

    M = datareg.toMNI;
    [U, L, V]  = svd(M(1:3, 1:3));
    M(1:3,1:3) = U*V';
    
    switch space
        case 'MNI-aligned'            
            BF.data.MEG.vol  = ft_transform_vol(M, vol);
            BF.data.MEG.sens = ft_transform_sens(M, sens);            
            
            BF.data.transforms.toMNI         = datareg.toMNI/M;
            BF.data.transforms.toMNI_aligned = eye(4);
            BF.data.transforms.toHead        = inv(M);
            BF.data.transforms.toNative      = D.inv{val}.mesh.Affine\BF.data.transforms.toMNI;
        case 'Head'
            BF.data.MEG.vol  = vol;
            BF.data.MEG.sens = sens;     
            
            BF.data.transforms.toMNI         = datareg.toMNI;
            BF.data.transforms.toMNI_aligned = M;
            BF.data.transforms.toHead        = eye(4);
            BF.data.transforms.toNative      = D.inv{val}.mesh.Affine\BF.data.transforms.toMNI;
        case 'Native'
            if istemplate
                error('Cannot work in Native space with template head, use MNI-aligned.');
            end
            
            BF.data.MEG.vol  = ft_transform_vol(D.inv{val}.mesh.Affine\datareg.toMNI, vol);
            BF.data.MEG.sens = ft_transform_vol(D.inv{val}.mesh.Affine\datareg.toMNI, sens);
            
            BF.data.transforms.toMNI         = D.inv{val}.mesh.Affine;
            BF.data.transforms.toMNI_aligned = M*datareg.fromMNI*BF.data.transforms.toMNI;
            BF.data.transforms.toHead        = datareg.fromMNI*BF.data.transforms.toMNI;
            BF.data.transforms.toNative      = eye(4);
    end
end
            

if eegind > 0
    vol      = D.inv{val}.forward(eegind).vol;
    datareg  = D.inv{val}.datareg(eegind);
    sens     = datareg.sensors;
    
    BF.data.EEG.vol  = vol;        
    BF.data.EEG.sens = sens;             
        
    if isfield(BF.data, 'transforms')  % With MEG
        if ~istemplate && isequal(space, 'Native')
            BF.data.EEG.vol  = vol;
            BF.data.EEG.sens = sens;
        elseif istemplate
            error('Combining EEG and MEG cannot be done with template head for now.');
        else
            if isa(vol, 'char')
                vol = ft_read_vol(vol);
            end
            
            BF.data.EEG.vol  = ft_transform_vol(inv(BF.data.transforms.toNative), vol);
            BF.data.EEG.sens = ft_transform_sens(inv(BF.data.transforms.toNative), sens);
        end
    else                             % EEG only        
        M = datareg.toMNI;
        [U, L, V]  = svd(M(1:3, 1:3));
        M(1:3,1:3) = U*V';
        
        
        switch space
            case 'Native'
                BF.data.EEG.vol  = vol;
                BF.data.EEG.sens = sens;
                
                BF.data.transforms.toMNI         = datareg.toMNI;
                BF.data.transforms.toMNI_aligned = M;
                BF.data.transforms.toHead        = eye(4); % Lets define Native and Head 
                BF.data.transforms.toNative      = eye(4); % to be the same thing in this case
            case {'MNI-aligned'}
                if isa(vol, 'char')
                    vol = ft_read_vol(vol);
                end
                
                BF.data.EEG.vol  = ft_transform_vol(M, vol);
                BF.data.EEG.sens = ft_transform_sens(M, sens);
                
                BF.data.transforms.toMNI         = datareg.toMNI/M;
                BF.data.transforms.toMNI_aligned = eye(4);
                BF.data.transforms.toHead        = inv(M);
                BF.data.transforms.toNative      = inv(M);
            case {'Head'}              
               error('Head space is not defined for EEG data');
        end
    end
end

BF.data.space = space;

bf_save(BF, 'overwrite');

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
