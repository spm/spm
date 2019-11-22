function res = bf_output_image_dics(BF, S)
% Computes DICS image
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_output_image_dics.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    all = cfg_const;
    all.tag = 'all';
    all.name = 'All';
    all.val  = {1};
    
    condlabel = cfg_entry;
    condlabel.tag = 'condlabel';
    condlabel.name = 'Condition label';
    condlabel.strtype = 's';
    condlabel.val = {''};
    
    conditions = cfg_repeat;
    conditions.tag = 'conditions';
    conditions.name = 'Conditions';
    conditions.help = {'Specify the labels of the conditions to be included in the inversion'};
    conditions.num  = [1 Inf];
    conditions.values  = {condlabel};
    conditions.val = {condlabel};
    
    whatconditions = cfg_choice;
    whatconditions.tag = 'whatconditions';
    whatconditions.name = 'What conditions to include?';
    whatconditions.values = {all, conditions};
    whatconditions.val = {all};
    
    sametrials = cfg_menu;
    sametrials.tag = 'sametrials';
    sametrials.name = 'Trials same as for filters';
    sametrials.labels = {'yes', 'no'};
    sametrials.values = {true, false};
    sametrials.val = {false};
    sametrials.help = {'Take the same trials as used for filter computation',...
        'This is useful for bootstrap.'};
    
    woi = cfg_entry;
    woi.tag = 'woi';
    woi.name = 'Time windows of interest';
    woi.strtype = 'r';
    woi.num = [Inf 2];
    woi.val = {[-Inf Inf]};
    woi.help = {'Time windows (in ms)'};
    
    name = cfg_entry;
    name.tag = 'name';
    name.name = 'Channel name';
    name.strtype = 's';
    name.num = [1 Inf];
    name.help = {'Reference channel name.'};
    
    shuffle         = cfg_menu;
    shuffle.tag     = 'shuffle';
    shuffle.name    = 'Shuffle';
    shuffle.help    = {'Shuffle the reference channel to produce the null case.'};
    shuffle.labels  = {'yes', 'no'};
    shuffle.values  = {1, 0};
    shuffle.val = {0};
    
    refchan      = cfg_branch;
    refchan.tag  = 'refchan';
    refchan.name = 'Reference channel';
    refchan.val  = {name, shuffle};
    
    refdip = cfg_entry;
    refdip.tag = 'refdip';
    refdip.name = 'Reference source';
    refdip.strtype = 'r';
    refdip.num = [1 3];
    refdip.help = {'Location of the reference in MNI coordinates'};
    
    power = cfg_const;
    power.tag = 'power';
    power.name = 'Power (no reference)';
    power.val  = {1};
    power.help = {'Compute power image'};
    
    reference = cfg_choice;
    reference.tag = 'reference';
    reference.name = 'Reference type';
    reference.values = {power, refchan, refdip};
    reference.val = {power};
    
    powmethod         = cfg_menu;
    powmethod.tag     = 'powmethod';
    powmethod.name    = 'Power summary method';
    powmethod.help    = {'How to summarise the power for vector beamformer'};
    powmethod.labels  = {'trace', 'lambda1', 'imag'};
    powmethod.values  = {'trace', 'lambda1', 'imag'};
    powmethod.val = {'lambda1'};
    
    foi = cfg_entry;
    foi.tag = 'foi';
    foi.name = 'Frequency band of interest';
    foi.strtype = 'r';
    foi.num = [1 2];
    foi.help = {'Frequency window within which to compute CSD over (Hz)'};
    
    taper = cfg_menu;
    taper.tag = 'taper';
    taper.name = 'Taper';
    taper.help = {'Save taper as well as power'};
    taper.labels = {'Hanning', 'Rectangular', 'DPSS', 'Sine'};
    taper.values = {'hanning', 'rectwin', 'dpss', 'sine'};
    taper.val = {'dpss'};
    
    contrast = cfg_entry;
    contrast.tag = 'contrast';
    contrast.name = 'Time contrast';
    contrast.strtype = 'r';
    contrast.num = [1 Inf];
    contrast.val = {1};
    
    logpower = cfg_menu;
    logpower.tag = 'logpower';
    logpower.name = 'Take log of power';
    logpower.labels = {'yes', 'no'};
    logpower.values = {true, false};
    logpower.val = {false};
    logpower.help = {'Take the log of power before computing time contrast',...
        'This is equivalent to log of the ratio.'};
    
    result         = cfg_menu;
    result.tag     = 'result';
    result.name    = 'What to output';
    result.help    = {'Specify output type.'};
    result.labels  = {
        'Single image'
        'Image per condition'
        'Image per trial'
        }';
    result.values  = {
        'singleimage'
        'bycondition'
        'bytrial'
        }';
    result.val = {'singleimage'};
    
    scale         = cfg_menu;
    scale.tag     = 'scale';
    scale.name    = 'Scale by filter norm';
    scale.help    = {'Scale by IID noise projected through the filters.'};
    scale.labels  = {'yes', 'no'};
    scale.values  = {'yes', 'no'};
    scale.val = {'no'};
    
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
    
    image_dics      = cfg_branch;
    image_dics.tag  = 'image_dics';
    image_dics.name = 'DICS image';
    image_dics.val  = {reference, powmethod, whatconditions, sametrials, woi,  contrast, logpower, foi, taper, result, scale, modality};
    
    res = image_dics;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

S.woi = 1e-3*S.woi; % ms -> s

samples = {};
for i = 1:size(S.woi, 1)
    samples{i} = D.indsample(S.woi(i, 1)):D.indsample(S.woi(i, 2));
end

if isfield(S.whatconditions, 'all')
    S.whatconditions.condlabel = D.condlist;
end

for i = 1:numel(S.whatconditions.condlabel)
    if S.sametrials
        trials{i} = BF.features.trials(strmatch(S.whatconditions.condlabel{i},...
            D.conditions(BF.features.trials)));
    else
        trials{i} = D.indtrial(S.whatconditions.condlabel{i}, 'GOOD');
    end
    
    if isempty(trials{i})
        error('No trials matched the selection.');
    end
    
end

if isempty(trials)
    error('No trials matched the selection, check the specified condition labels');
end

alltrials = spm_vec(trials);
ntrials   = length(alltrials);

channels = BF.features.(S.modality).chanind;
U        = BF.features.(S.modality).U;
nchan    = size(U, 1);

Cf  = {};
refindx = [];
Wr = [];
if isfield(S.reference, 'refchan')
    refindx = D.indchannel(S.reference.refchan.name);
    
    if isempty(refindx)
        error('The specified reference channel was not found in the data.');
    end
    
    Cr = {};
    Pr = [];
    
    if S.reference.refchan.shuffle
        shuffle = randperm(ntrials);
        prefix = 'dics_refcoh_shuffled';
    else
        shuffle = 1:ntrials;
        prefix = 'dics_refcoh';
    end
elseif isfield(S.reference, 'refdip')
    % transform coords in MNI space into space where we are doing the beamforming
    seed = spm_eeg_inv_transform_points(inv(BF.data.transforms.toMNI), S.reference.refdip);
    pos  = BF.sources.pos;
    nvert = size(pos, 1);
    
    dist = sqrt(sum((pos - repmat(seed, nvert, 1)).^2, 2));
    
    [mdist, ind] = min(dist);
    
    if mdist > 20
        warning(['Closest match is ' mdist ' mm away from the specified location.']);
    end
    
    Wr =  BF.inverse.(S.modality).W{ind};
    
    prefix = 'dics_dipcoh';
else
    prefix = 'dics_pow';
end


nsamples = unique(cellfun(@length, samples));
if length(nsamples) > 1
    error('All time windows should be equal lentgh')
end

centerfreq = mean(S.foi);
tapsmofrq  = 0.5*(abs(diff(S.foi)));

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials , 'Computing CSD'); drawnow;
if ntrials  > 100, Ibar = floor(linspace(1, ntrials ,100));
else Ibar = 1:ntrials; end

for i = 1:ntrials
    for j = 1:numel(samples)
        Y  = U'*squeeze(D(channels, samples{j}, alltrials(i)));
        
        [fourier, ntap] = ft_specest_mtmfft(Y, D.time(samples{j}), 'freqoi', centerfreq, ...
            'tapsmofrq', tapsmofrq, 'taper', S.taper, 'verbose', 0);
        
        dat  = transpose(fourier);
        
        Cf{i, j} = (dat * ctranspose(dat)) ./ ntap;
        
        if ~isempty(refindx)
            Y  = squeeze(D(refindx, samples{j}, alltrials(shuffle(i))));
            
            [fourier, ntap] = ft_specest_mtmfft(Y, D.time(samples{j}), 'freqoi', centerfreq, ...
                'tapsmofrq', tapsmofrq, 'taper', S.taper, 'verbose', 0);
            
            ref = transpose(fourier);
            Cr{i, j}  = dat * ctranspose(ref) ./ ntap;
            Pr(i, j)  = ref * ctranspose(ref) ./ ntap;
        end
    end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

W = BF.inverse.(S.modality).W;
nvert = numel(W);

mCf = {};

if ~isempty(refindx)
    mCr = {};
    mPr = [];
end

condind = spm_unvec(1:ntrials, trials);

switch S.result
    case 'singleimage'
        for i = 1:size(Cf, 2)
            mCf{1, i} = squeeze(sum(cat(3, Cf{:, i}), 3))./ntrials;
            
            if ~isempty(refindx)
                mCr{1, i} = squeeze(sum(cat(2, Cr{:, i}), 2))./ntrials;
                mPr(1, i) = sum(Pr(:, i))./ntrials;
            end
        end
    case 'bycondition';
        for c = 1:numel(condind)
            for i = 1:size(Cf, 2)
                mCf{c, i} = squeeze(sum(cat(3, Cf{condind{c}, i}), 3))./length(condind{c});
                
                if ~isempty(refindx)
                    mCr{c, i} = squeeze(sum(cat(2, Cr{condind{c}, i}), 2))./length(condind{c});
                    mPr(c, i) = sum(Pr(condind{c}, i))./length(condind{c});
                end
            end
        end
    case 'bytrial'
        for c = 1:ntrials
            mCf = Cf;
            
            if ~isempty(refindx)
                error('Coherence image cannot be computed for single trials');
            end
            
        end
end

spm('Pointer', 'Watch');drawnow;

for c = 1:size(mCf, 1)
    spm_progress_bar('Init', nvert, ...
        sprintf('Scanning grid points image %d/%d', c, size(mCf, 1))); drawnow;
    if nvert > 100, Ibar = floor(linspace(1, nvert,100));
    else Ibar = 1:nvert; end
    
    pow = nan(1, nvert);
    cpow = zeros(1, numel(mCf(c, :)));
    
    for i = 1:nvert
        if ~isnan(W{i})
            w    = W{i};
            
            for j = 1:numel(mCf(c, :))
                if ~isempty(refindx)
                    estimate = ft_inverse_beamformer_dics(w, mCf{c, j}, 'Cr', mCr{c, j}, 'Pr', mPr(c, j), ...
                        'filterinput', 'yes',  'projectnoise', S.scale, 'powmethod', S.powmethod, ...
                        'keepfilter', 'no', 'keepleadfield', 'no', 'keepcsd', 'no', 'feedback', 'none');
                    
                    cpow(j) = estimate.coh;
                elseif ~isempty(Wr)
                    estimate = ft_inverse_beamformer_dics(w, mCf{c, j}, 'refdip', Wr, ...
                        'filterinput', 'yes',  'projectnoise', S.scale, 'powmethod', S.powmethod, ...
                        'keepfilter', 'no', 'keepleadfield', 'no', 'keepcsd', 'no', 'feedback', 'none');
                    
                    cpow(j) = estimate.coh;
                else
                    estimate = ft_inverse_beamformer_dics(w, mCf{c, j},...
                        'filterinput', 'yes',  'projectnoise', S.scale, 'powmethod', S.powmethod, ...
                        'keepfilter', 'no', 'keepleadfield', 'no', 'keepcsd', 'no', 'feedback', 'none');
                    
                    cpow(j) = estimate.pow;
                end
                
                if isequal(S.scale, 'yes')
                    cpow(j) = cpow(j)./estimate.noise;
                end
                
                if S.logpower
                    pow(i) = log(cpow)*S.contrast';
                else
                    pow(i) = cpow*S.contrast';
                end
                
            end
        end
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end
    end
    
    spm_progress_bar('Clear');
    
    image(c).val   = pow;
    
    switch S.result
        case 'singleimage'
            image(c).label = [prefix '_' spm_file(D.fname, 'basename')];
        case 'bycondition'
            image(c).label = [prefix '_cond_' S.whatconditions.condlabel{c} '_' spm_file(D.fname, 'basename')];
        case 'bytrial'
            for k = 1:numel(condind)
                if any(c == condind{k})
                    break;
                end
            end
            image(c).label = [prefix '_cond_' S.whatconditions.condlabel{k}...
                '_trial_' num2str(alltrials(c)) '_' spm_file(D.fname, 'basename')];
    end
end

spm('Pointer', 'Arrow');drawnow;

res = image;
