function res = bf_inverse_deflect(BF, S)
% Used DeFleCT framework to compute spatial filters.
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, using the code from Matti Stenroos and Olaf Hauk
% http://imaging.mrc-cbu.cam.ac.uk/meg/AnalyzingData/DeFleCT_SpatialFiltering_Tools
%
% Please cite:
% Hauk O, Stenroos M.
% A framework for the design of flexible cross-talk functions for spatial filtering of EEG/MEG data: DeFleCT.
% Human Brain Mapping 2013
% $Id: bf_inverse_deflect.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    
    list = cfg_entry;
    list.tag = 'list';
    list.name = 'List vertex indices';
    list.strtype = 'n';
    list.num = [1 Inf];
    list.help = {'Specify sources of interest by listing vertex indices.'};
    
    pos = cfg_entry;
    pos.tag = 'pos';
    pos.name = 'MNI coordinates';
    pos.strtype = 'r';
    pos.num = [1 3];
    pos.help = {'Locations for the VOI center in MNI coordinates'};
    pos.val = {};
    
    radius = cfg_entry;
    radius.tag = 'radius';
    radius.name = 'Radius';
    radius.strtype = 'r';
    radius.num = [1 1];
    radius.val = {0};
    radius.help = {'Radius (in mm) for the VOI (leave 0 for closest point)'};
    
    voi = cfg_branch;
    voi.tag = 'voi';
    voi.name = 'VOI';
    voi.val = {pos, radius};
    
    passband = cfg_repeat;
    passband.tag = 'passband';
    passband.name = 'Passband sources';
    passband.num  = [1 Inf];
    passband.values = {voi, list};
    
    stopband = cfg_repeat;
    stopband.tag = 'stopband';
    stopband.name = 'Stopband sources';
    stopband.num  = [0 Inf];
    stopband.values = {voi, list};
    
    svdpassband = cfg_entry;
    svdpassband.tag = 'svdpassband';
    svdpassband.name = 'SVD passband';
    svdpassband.strtype = 'w';
    svdpassband.num = [1 1];
    svdpassband.val = {0};
    svdpassband.help = {'Number of components to summarise the passband in.,',...
        'Leave at zero for no SVD.'};
    
    svdstopband = cfg_entry;
    svdstopband.tag = 'svdstopband';
    svdstopband.name = 'SVD stopband';
    svdstopband.strtype = 'w';
    svdstopband.num = [1 1];
    svdstopband.val = {0};
    svdstopband.help = {'Number of components to summarise the stopband in.,',...
        'Leave at zero for no SVD.'};
    
    forcepassband         = cfg_menu;
    forcepassband.tag     = 'forcepassband';
    forcepassband.name    = 'Force passband';
    forcepassband.help    = {'Forces the output for all passband components '};
    forcepassband.labels  = {'yes', 'no'};
    forcepassband.values  = {1, 0};
    forcepassband.val = {0};
    
    label = cfg_entry;
    label.tag = 'label';
    label.name = 'Label';
    label.strtype = 's';
    label.help = {'Label for the output source'};
    
    usecov         = cfg_menu;
    usecov.tag     = 'usecov';
    usecov.name    = 'Use covariance matrix';
    usecov.help    = {'Use covariance matrix for pre-whitening'};
    usecov.labels  = {'yes', 'no'};
    usecov.values  = {1, 0};
    usecov.val = {1};
    
    filter = cfg_branch;
    filter.tag = 'filter';
    filter.name = 'Filter';
    filter.val = {label, passband,  svdpassband, forcepassband, stopband, svdstopband, usecov};
    
    filters = cfg_repeat;
    filters.tag = 'filters';
    filters.name = 'Filters';
    filters.num  = [1 Inf];
    filters.values = {filter};    
    
    snr = cfg_entry;
    snr.tag = 'snr';
    snr.name = 'SNR';
    snr.strtype = 'r';
    snr.num = [1 1];
    snr.val = {5};
    snr.help = {'The assumed ratio of variances of signal and noise,',...
        'used for setting the regularisation parameter.'};
    
    trunc = cfg_entry;
    trunc.tag = 'trunc';
    trunc.name = 'Truncation parameter';
    trunc.strtype = 'w';
    trunc.num = [1 1];
    trunc.val = {0};
    trunc.help = {'The number of (smallest) singular values of the covariance matrix that are set to ',...
        'zero before making the whitener. For example, if the data has been SSP-projected, it needs to be at least the number of ',...
        'components projected away.'};
     
    
    deflect     = cfg_branch;
    deflect.tag  = 'deflect';
    deflect.name = 'DeFleCT';
    deflect.val  = {filters, snr, trunc};
    deflect.help = {'DeFleCT spatial filter design framework by Matti Stenroos and Olaf Hauk',...
        'http://imaging.mrc-cbu.cam.ac.uk/meg/AnalyzingData/DeFleCT_SpatialFiltering_Tools',...
        'Please cite:'...
        'Hauk O, Stenroos M.',...
        'A framework for the design of flexible cross-talk functions for spatial filtering of EEG/MEG data: DeFleCT.',...
        'Human Brain Mapping 2013'};
    res = deflect;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end


C     =  BF.features.(S.modality).C;
U     =  BF.features.(S.modality).U;

L = [];
Li = {};
for i = 1:numel(S.L)
    cL    = U'*S.L{i};
    Li{i} = size(L, 2)+(1:size(cL, 2));
    L  =[L cL];
end

mnipos = spm_eeg_inv_transform_points(BF.data.transforms.toMNI,  BF.sources.pos);

nfilt = numel(S.filter);

W     = cell(1, nfilt);
label = {};

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nfilt, ['Computing ' S.modality ' filters']); drawnow;
if nfilt > 100, Ibar = floor(linspace(1, nfilt,100));
else Ibar = 1:nfilt; end

for i = 1:nfilt
    filter = S.filter(i);
    if ~filter.svdpassband
        filter.svdpassband = [];
    end
    
     if ~filter.svdstopband
        filter.svdstopband = [];
    end
    
    passband = get_vertices(filter.passband, mnipos);
    stopband = get_vertices(filter.stopband, mnipos);
    
    passband = cat(2, Li{passband});
    
    stopband = cat(2, Li{stopband});
    
    if filter.usecov
        W{i} = DeFleCT(passband, filter.svdpassband, filter.forcepassband, stopband, filter.svdstopband,...
            L, C, S.snr, S.trunc);
    else
        W{i} = DeFleCT(passband, filter.svdpassband, filter.forcepassband, stopband, filter.svdstopband,...
            L, [], S.snr, eye(size(L, 1)));
    end    
    
    label{i} = filter.label;
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end


spm_progress_bar('Clear');

disp('Note: it is crucial to visualise the filter-cosstalk functions to verify correct design.');

res.W     = W;
res.label = label;


function vertices = get_vertices(band, mnipos)

vertices = [];
for j = 1:numel(band)
    switch char(fieldnames(band{j}))
        case 'voi'
            pnt = band{j}.voi.pos;
            rad = band{j}.voi.radius;
            
            dist = sqrt(sum((mnipos - repmat(pnt, size(mnipos, 1), 1)).^2, 2));
            
            if ~rad
                [mindist, ind] = min(dist);
                if mindist>30
                    error('There are no sources within 3cm of the specified location');
                end
            else
                ind = find(dist<=rad);
            end
        case 'list'
            ind = band{j}.list;
    end
    
    vertices = [vertices; ind(:)];
end