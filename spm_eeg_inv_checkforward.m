function spm_eeg_inv_checkforward(varargin)
% Check M/EEG forward model
% FORMAT spm_eeg_inv_checkforward(D, val, ind)
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_checkforward.m 5592 2013-07-24 16:25:55Z vladimir $


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
    sens     = D.inv{val}.datareg(ind).sensors;    
    Mcortex  = forward(ind).mesh;
catch
    warndlg('please coregister and create forward model')
    return
end

Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf

if ismac
    set(Fgraph,'renderer','zbuffer');
else
    set(Fgraph,'renderer','OpenGL');
end


spm('Pointer', 'Watch');drawnow;
%--------------------------------------------------------------------------
chanind = strmatch(modality, D.chantype);
chanind = setdiff(chanind, D.badchannels);
if isempty(chanind)
    error(['No good ' modality ' channels were found.']);
end

if ischar(vol)
    vol = ft_read_vol(vol);
end

face    = Mcortex.face;
vert    = Mcortex.vert;
h_ctx   = patch('vertices',vert,'faces',face,'EdgeColor','b','FaceColor','b');

hold on

ft_plot_vol(vol, 'edgecolor', [0 0 0], 'facealpha', 0);

ft_plot_sens(sens, 'style', '*b');

rotate3d on;

axis off
axis vis3d
axis equal

spm('Pointer', 'Arrow');
