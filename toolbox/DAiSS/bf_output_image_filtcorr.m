function res = bf_output_image_filtcorr(BF, S)
% Computes filter correlation images
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_output_image_filtcorr.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    pos = cfg_entry;
    pos.tag = 'pos';
    pos.name = 'Seed MNI coordinates';
    pos.strtype = 'r';
    pos.num = [1 3];
    pos.help = {'Locations for the seed in MNI coordinates (closest point is chosen'};
    pos.val = {};
    
    label = cfg_entry;
    label.tag = 'label';
    label.name = 'Label';
    label.strtype = 's';
    label.help = {'Label for source of interest'};    
    
    seedspec = cfg_choice;
    seedspec.tag = 'seedspec';
    seedspec.name = 'Seed specification';
    seedspec.values = {pos, label};
    
    corrtype         = cfg_menu;
    corrtype.tag     = 'corrtype';
    corrtype.name    = 'Correlation type';
    corrtype.help    = {'Whether to correlate filters with other filters of with other leadfields'};
    corrtype.labels  = {'Filter-Filter', 'Filter-Leadfield'};
    corrtype.values  = {'filtfilt', 'filtlf'};
    corrtype.val = {'filtfilt'};
    
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
    
    filtcorr      = cfg_branch;
    filtcorr.tag  = 'image_filtcorr';
    filtcorr.name = 'Filter correlations image';
    filtcorr.val  = {seedspec, corrtype, modality};
    
    res = filtcorr;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

nvert = size(BF.sources.pos, 1);

% transform coords in MNI space into space where we are doing the beamforming
if isfield(S.seedspec, 'pos')
    mnipos  = spm_eeg_inv_transform_points(BF.data.transforms.toMNI, BF.sources.pos);  
    
    dist = sqrt(sum((mnipos - repmat(S.seedspec.pos, nvert, 1)).^2, 2));
    
    [mdist, ind] = min(dist);
    
    if mdist > 20
        warning(['Closest match is ' mdist ' mm away from the specified location.']);
    end    
else
    if isfield(BF.inverse.(S.modality), 'label')
        ind = strmatch(S.seedspec.label, BF.inverse.(S.modality).label, 'exact');
    else
        error('Filters are not labeled, use position to specify seed.');
    end    
end

ws =  BF.inverse.(S.modality).W{ind};

[QA, dum] = qr(orth(ws'),0);
Q         = svd(QA'*QA);
scale     = 1/sum(Q);

spm('Pointer', 'Watch');drawnow;

spm_progress_bar('Init', nvert, 'Scanning grid points'); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

pow = nan(1, nvert);

U     =  BF.features.(S.modality).U;

for i = 1:nvert
    switch S.corrtype
        case 'filtfilt'
            w = BF.inverse.(S.modality).W{i}';
        case 'filtlf'
            w =  U'*BF.inverse.(S.modality).L{i};
    end
    
    if ~isnan(w)
        % This is subspace intersection which is supposed to handle the case
        % when filters are more than 1D and be equivalent to correlation
        % coefficient for the 1D case        
        [QB, dum] = qr(orth(w),0);
        Q = svd(QA'*QB);
        pow(i) = scale*sum(Q);
    end
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

image.val   = pow;

image.label = ['filtcorr_'  spm_file(fname(BF.data.D), 'basename')];
 
spm('Pointer', 'Arrow');drawnow;

res = image;