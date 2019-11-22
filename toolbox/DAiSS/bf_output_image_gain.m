function res = bf_output_image_gain(BF, S)
% Computes gain image
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Ashwini Oswal, Vladimir Litvak
% $Id: bf_output_image_gain.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0  
    type         = cfg_menu;
    type.tag     = 'type';
    type.name    = 'Output type';
    type.help    = {'The calculation that should be performed'};
    type.labels  = {'Reduced/Original'};
    type.values  = {'reduced_vs_orig'};
    type.val     = {'reduced_vs_orig'};
    
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
    
    gain      = cfg_branch;
    gain.tag  = 'image_gain';
    gain.name = 'Gain image';
    gain.val  = {type, modality};
    
    res = gain;
    
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

U     =  BF.features.(S.modality).U;

for i = 1:nvert
    switch S.type
        case 'reduced_vs_orig'
            lt = U;
            
            if isfield(BF, 'inverse')
                lf = BF.inverse.(S.modality).L{i};
            else
                lf = BF.sources.L.(S.modality(1:3)){i};
            end
            if ~isnan(lf)
                % spatial filtering function defined in equation (7) of VanVeen
                % paper NB there is a slight error in the paper itself
                val(i) = trace(lf'*lt*(lf'*lt)')/trace(lf'*lf);
            end
    end    
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

image.val   = val;

image.label = ['gain_'  spm_file(fname(BF.data.D), 'basename')];
 
spm('Pointer', 'Arrow');drawnow;

res = image;