function res = bf_write_gifti(BF, S)
% Writes out beamformer results as GIfTI meshes
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_write_gifti.m 7703 2019-11-22 12:06:29Z guillaume $

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
    
    visualise = cfg_menu;
    visualise.tag = 'visualise';
    visualise.name = 'Visualise the outputs';
    visualise.help = {'Visualise the outputs in the graphics window.'};
    visualise.labels = {'inflated', 'yes', 'no'};
    visualise.values = {2, 1, 0};
    visualise.val = {0};
    
    res      = cfg_branch;
    res.tag  = 'gifti';
    res.name = 'GIfTI';
    res.val  = {normalise, space, visualise};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

if ~isfield(BF.sources, 'mesh')
    error('Only mesh source space is supported');
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
        source = BF.sources.mesh.canonical;
    case 'aligned'
        source = BF.sources.mesh.individual;
        source.vert =  spm_eeg_inv_transform_points(BF.data.transforms.toMNI_aligned, source.vert);
    case 'native'
        source = BF.sources.mesh.individual;
        source.vert =  spm_eeg_inv_transform_points(BF.data.transforms.toNative, source.vert);
end

save(gifti(source), [BF.data.D.fname '.surf.gii']);

nimages = numel(BF.output.image);

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nimages , 'Writing out images'); drawnow;
if nimages  > 100, Ibar = floor(linspace(1, nimages ,100));
else Ibar = 1:nimages; end

for i = 1:nimages
    fname = fullfile(pwd, [BF.output.image(i).label '.gii']);
    
    G = gifti;
    G.private.metadata(1).name = 'SurfaceID';
    G.private.metadata(1).value = [BF.data.D.fname '.surf.gii'];
    
    G.cdata = scale(i)*BF.output.image(i).val;
    
    G.cdata = G.cdata(:);
    
    save(G, fname, 'ExternalFileBinary');
            
    res.files{i, 1} = fname;
    
    if S.visualise
        Fgraph  = spm_figure('GetWin', fname); figure(Fgraph); clf
        a = axes;
        if S.visualise == 2
            H = spm_mesh_render('Disp', spm_mesh_inflate(gifti([BF.data.D.fname '.surf.gii'])), 'Parent', a);
            H = spm_mesh_render('Overlay', H, gifti(fname));
        else
            H = spm_mesh_render('Disp',gifti([BF.data.D.fname '.surf.gii']), 'Parent', a);
            H = spm_mesh_render('Overlay', H, gifti(fname));
        end
         spm_mesh_render('Colormap', H, jet(256));
    end
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');
