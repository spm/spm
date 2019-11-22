function spm_eeg_inv_checkforward(varargin)
% Check M/EEG forward model
% FORMAT spm_eeg_inv_checkforward(D, val, ind)
%__________________________________________________________________________
% Copyright (C) 2008-2018 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_checkforward.m 7702 2019-11-22 11:32:26Z guillaume $


%-SPM data structure
%--------------------------------------------------------------------------
[D, val] = spm_eeg_inv_check(varargin{:});

forward = D.inv{val}.forward;

if nargin < 3
    str = sprintf('%s|', forward(:).modality);
    str = str(1:(end-1));
    
    ind = spm_input('What to display?','+1','b',str,1:numel(forward),1);
else
    ind = varargin{3};
end

try
    vol      = forward(ind).vol;
    modality = forward(ind).modality;
    if isfield(forward(ind), 'siunits') && forward(ind).siunits
        sens     = forward(ind).sensors;
    else %backward compatibility
        sens     = D.inv{val}.datareg(ind).sensors;
    end
    Mcortex  = forward(ind).mesh;
catch
    spm('alert!','Please coregister and create forward model',mfilename);
    return
end

%-Display
%--------------------------------------------------------------------------
Fgraph  = spm_figure('GetWin','Graphics');
spm_figure('Focus',Fgraph);
spm_figure('Clear',Fgraph);

spm('Pointer', 'Watch');

chanind = strmatch(modality, D.chantype);
chanind = setdiff(chanind, D.badchannels);
if isempty(chanind)
    error(['No good ' modality ' channels were found.']);
end

if ischar(vol)
    vol = ft_read_headmodel(vol);
end

face    = Mcortex.face;
vert    = Mcortex.vert;
h_ctx   = patch('vertices',vert,'faces',face,'EdgeColor','b','FaceColor','b');

hold on

[volp, sens] = ft_prepare_vol_sens(vol, sens, 'channel', D.chanlabels(chanind));

ft_plot_headmodel(vol, 'edgecolor', [0 0 0], 'facealpha', 0);

hold on

try
    ft_plot_sens(sens, 'style', '*', 'edgecolor', 'g', 'elecsize', 20, 'coil', ft_senstype(sens, 'eeg'));
catch
    ft_plot_sens(sens, 'edgecolor', 'g', 'coilshape', 'point', 'elecsize', 20, 'coil', true);
end

rotate3d on;

axis off
axis vis3d
axis equal

spm('Pointer', 'Arrow');
