function res = bf_write_nifti(BF, S)
% Writes out nifti images of beamformer results
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_write_nifti.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    normalise         = cfg_menu;
    normalise.tag     = 'normalise';
    normalise.name    = 'Global normalisation';
    normalise.help    = {'Normalise image values by the mean'};
    normalise.labels  = {
        'No'
        'Each image separately'
        'Across images'
        }';
    normalise.values  = {
        'no'
        'separate'
        'all'
        }';
    normalise.val = {'no'};
    
    space         = cfg_menu;
    space.tag     = 'space';
    space.name    = 'Image space';
    space.help    = {'Specify image space'};
    space.labels  = {
        'MNI'
        'Native'
        'MNI-aligned'
        }';
    space.values  = {
        'mni'
        'native'
        'aligned'
        }';
    space.val = {'mni'};
    
    nifti      = cfg_branch;
    nifti.tag  = 'nifti';
    nifti.name = 'NIfTI';
    nifti.val  = {normalise, space};
    
    res = nifti;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

scale = ones(1, numel(BF.output.image));
switch S.normalise
    case  'separate'
        for i = 1:numel(BF.output.image)
            val = BF.output.image(i).val;
            scale(i) = 1./mean(abs(val(~isnan(val))));
        end
    case  'all'
        val = spm_vec({BF.output.image(:).val});
        scale = scale./mean(abs(val(~isnan(val))));
end

switch S.space
    case 'mni'
        sMRI   = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
    case 'aligned'
        sMRI   = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
    case 'native'
        sMRI   = BF.data.mesh.sMRI;
end

[pth,nam,ext,num] = spm_fileparts(sMRI);
sMRI = fullfile(pth, [nam ext]);

if isfield(BF.sources, 'grid') || isfield(BF.sources, 'voi')       
    if isfield(BF.sources, 'grid')
        sourcespace    = 'grid';
        source         = BF.sources.grid;
        source.pos     = BF.sources.grid.allpos;
    else
        sourcespace    = 'voi';
        source         = [];
        source.pos     = BF.sources.pos;
        source.inside  = 1:size(BF.sources.pos, 1);
        source.outside = [];
    end
    

    
    switch S.space
        case 'mni'
            source = ft_transform_geometry(BF.data.transforms.toMNI, source);
        case 'aligned'
            source = ft_transform_geometry(BF.data.transforms.toMNI_aligned, source);
        case 'native'
            source = ft_transform_geometry(BF.data.transforms.toNative, source);
    end
    
    cfg = [];
    cfg.parameter = 'pow';
    cfg.downsample = 1;
    cfg.showcallinfo = 'no';
    
elseif isfield(BF.sources, 'mesh')
    sourcespace = 'mesh';
    
    switch S.space
        case 'mni'
            source = BF.sources.mesh.canonical;
        case 'aligned'
            source = BF.sources.mesh.individual;
            source.vert =  spm_eeg_inv_transform_points(BF.data.transforms.toMNI_aligned, source.vert);
        case 'native'
            source = BF.sources.mesh.individual;
            source.vert =  spm_eeg_inv_transform_points(BF.data.transforms.toNative, source.vert);
    end
    
    source = export(gifti(source), 'patch');
else
    error('Unsupported source space type');
end


outvol = spm_vol(sMRI);
outvol.dt(1) = spm_type('float32');


nimages = numel(BF.output.image);

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nimages , 'Writing out images'); drawnow;
if nimages  > 100, Ibar = floor(linspace(1, nimages ,100));
else Ibar = 1:nimages; end


for i = 1:nimages
    outvol.fname= fullfile(pwd, [BF.output.image(i).label '.nii']);
    outvol = spm_create_vol(outvol);
    
    source.pow = scale(i)*BF.output.image(i).val;
    
    source.pow = source.pow(:);
    
    switch sourcespace
        case 'grid'
            pow = source.pow;
            source.pow = nan(size(source.pos, 1), 1);
            source.pow(source.inside) = pow;
            sourceint = ft_sourceinterpolate(cfg, source, ft_read_mri(sMRI, 'dataformat', 'nifti_spm'));
            Y = sourceint.pow;
        case 'mesh'
            Y = spm_mesh_to_grid(source, outvol, source.pow);
            spm_smooth(Y, Y, 1);
            Y = Y.*(Y > max(source.pow)*exp(-8));
        case 'voi'
            cfg.interpmethod = 'sphere_avg';
            cfg.sphereradius = 5;
            sourceint = ft_sourceinterpolate(cfg, source, ft_read_mri(sMRI, 'dataformat', 'nifti_spm'));
            Y = sourceint.pow;
            Y = reshape(Y, sourceint.dim);
    end
    
    spm_write_vol(outvol, Y);
    
    nifti.files{i, 1} = outvol.fname;
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

res = nifti;
