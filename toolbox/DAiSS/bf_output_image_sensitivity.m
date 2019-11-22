function res = bf_output_image_sensitivity(BF, S)
% Sensitivity profile for a group of sensors
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_output_image_sensitivity.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0      
    modality         = cfg_menu;
    modality.tag     = 'modality';
    modality.name    = 'Modality';
    modality.help    = {'Specify modality'};
    modality.labels  = {
        'MEG'
        'MEGPLANAR'
        'EEG'
        }';
    modality.values  = {
        'MEG'
        'MEGPLANAR'
        'EEG'
        }';
    modality.val = {'MEG'};
    
    sensitivity      = cfg_branch;
    sensitivity.tag  = 'image_sensitivity';
    sensitivity.name = 'Sensitivity image';
    sensitivity.val  = {spm_cfg_eeg_channel_selector, modality};
    
    res = sensitivity;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

nvert = size(BF.sources.pos, 1);

spm('Pointer', 'Watch');drawnow;

spm_progress_bar('Init', nvert, 'Scanning grid points'); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

val = nan(1, nvert);
D = BF.data.D;

selectchan = D.chanlabels(D.selectchannels(spm_cfg_eeg_channel_selector(S.channels)));

for i = 1:nvert   
    
    if isfield(BF, 'inverse')
        lf = BF.inverse.(S.modality).L{i};
        channels = BF.inverse.(S.modality).channels;	
    else
        lf = BF.sources.L.(S.modality(1:3)){i};
        channels = BF.sources.channels.(S.modality(1:3));
    end
    
    [sel1, sel2] = spm_match_str(selectchan, channels);
    
    lf = lf(sel2, :);
    
    if ~any(isnan(lf))        
        val(i) = sqrt(trace(lf'*lf));
    end
    
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

image.val   = val;

image.label = ['sensitivity_'  spm_file(fname(BF.data.D), 'basename')];
 
spm('Pointer', 'Arrow');drawnow;

res = image;