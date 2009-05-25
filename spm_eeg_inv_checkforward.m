function spm_eeg_inv_checkforward(varargin)
% checks forward model
% FORMAT spm_eeg_inv_checkforward(D)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_checkforward.m 3145 2009-05-25 16:00:43Z vladimir $

% SPM data structure
%==========================================================================
[D, val] = spm_eeg_inv_check(varargin{:});

forward = D.inv{val}.forward;

str = sprintf('%s|', forward(:).modality);
str = str(1:(end-1));

ind = spm_input('What to display?','+1', 'b',  str, 1:numel(forward), 1);    

try
    vol = forward(ind).vol;
    modality = forward(ind).modality;
    sens = D.inv{val}.datareg(ind).sensors;    
    Mcortex = forward(ind).mesh;
catch
    warndlg('please coregister and create forward model')
    return
end

Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf
spm('Pointer', 'Watch');drawnow;
%--------------------------------------------------------------------------
chanind = strmatch(modality, D.chantype);
chanind = setdiff(chanind, D.badchannels);
if isempty(chanind)
    error(['No good ' modality ' channels were found.']);
end

if ischar(vol)
    vol = fileio_read_vol(vol);
end

%--------------------------------------------------------------------------
cfg = [];

switch modality
    case 'EEG'
        cfg.elec = sens;
    case 'MEG'
        cfg.grad = sens;
    otherwise
        error('Unsupported modality');
end

cfg.channel           = D.chanlabels(chanind);
cfg.vol               = vol;
cfg.inwardshift       = 0;
cfg.plotgrid          = 'no';    
cfg.plotheadsurface   = 'no';
cfg.plotspheres       = 'yes';
cfg.plotbnd           = 'yes';
cfg.plotspherecenter  = 'yes'; 
cfg.plotlines         = 'no'; 
cfg.surftype          = 'edges';
cfg.surface_edgecolor = [0 0 0];

cfg.spheremesh        = 162;

ft_headmodelplot(cfg);

hold on

face    = Mcortex.face;
vert    = Mcortex.vert;
h_ctx   = patch('vertices',vert,'faces',face,'EdgeColor','b','FaceColor','b');

rotate3d on;
spm('Pointer', 'Arrow');drawnow;




