function res = bf_output_image_power(BF, S)
% Computes power image
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_output_image_power.m 7703 2019-11-22 12:06:29Z guillaume $

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
    woi.help = {'Time windows (in ms). N.B. Make sure these windows lie within the covariance window'};
    
    foi = cfg_entry;
    foi.tag = 'foi';
    foi.name = 'Frequency bands of interest';
    foi.strtype = 'r';
    foi.num = [Inf 2];
    foi.help = {'Frequency windows within which to power changes over (Hz) .N.B. Check this lies within covariance window'};
    
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
    scale.name    = 'Scale power by filter norm';
    scale.help    = {'Scale power by norm of the filters'};
    scale.labels  = {'norm only', 'with noise', 'no scaling'};
    scale.values  = {2, 1, 0};
    scale.val = {1};
    
    powermethod         = cfg_menu;
    powermethod.tag     = 'powermethod';
    powermethod.name    = 'Power summary method';
    powermethod.help    = {'How to summarise the power for vector beamformer'};
    powermethod.labels  = {'trace', 'lambda1'};
    powermethod.values  = {'trace', 'lambda1'};
    powermethod.val = {'trace'};
    
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
    
    image_power      = cfg_branch;
    image_power.tag  = 'image_power';
    image_power.name = 'Power image';
    image_power.val  = {whatconditions, sametrials, woi, foi, contrast, logpower, result, scale, powermethod, modality};
    
    res = image_power;
    
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

nsamples = unique(cellfun(@length, samples));

if length(nsamples)~=1,
    error('all windows must be of equal length');
end

windowduration  = nsamples/D.fsample;

dctfreq         = (0:nsamples-1)/2/windowduration;   % DCT frequencies (Hz)
dctT            = spm_dctmtx(nsamples,nsamples);

nbands = size(S.foi, 1);
allfreqind=[];

for fband = 1:nbands, %% allows one to break up spectrum and ignore some frequencies
    
    freqrange  = S.foi(fband,:);
    
    j          = find( (dctfreq >= freqrange(1)) & (dctfreq<=freqrange(2)));
  
    allfreqind = [allfreqind j];
    
end % for fband=1:Nbands

allfreqind = sort(unique(allfreqind));
if isempty(allfreqind),
    error('No valid frequency range found');
end;

Tband = dctT(:, allfreqind); % filter to this band

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

channels = BF.features.(S.modality).chanind;
U        = BF.features.(S.modality).U;
nchan    = size(U, 2);

YY       = {};
nsamples = unique(cellfun(@length, samples));
if length(nsamples) > 1
    error('All time windows should be equal lentgh')
end

alltrials = spm_vec(trials);
ntrials   = length(alltrials);

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials , 'Computing covariance'); drawnow;
if ntrials  > 100, Ibar = floor(linspace(1, ntrials ,100));
else Ibar = 1:ntrials; end

sumYY = 0;
N     = 0;
for i = 1:ntrials
    for j = 1:numel(samples)
        Y  = U'*squeeze(D(channels, samples{j}, alltrials(i)));
        Y  = detrend(Y', 'constant')';
        
        Y = Y*Tband;
        
        YY{i, j} = Y*Y';
        sumYY = sumYY+YY{i,j};
        N  = N+length(samples{j});
    end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end
sumYY = sumYY./N;

spm_progress_bar('Clear');

W = BF.inverse.(S.modality).W;
nvert = numel(W);

Cy = {};
condind = spm_unvec(1:ntrials, trials);

switch S.result
    case 'singleimage'
        for i = 1:size(YY, 2)
            Cy{1, i} = squeeze(sum(cat(3, YY{:, i}), 3))./((nsamples-1)*ntrials);
        end
    case 'bycondition';
        for c = 1:numel(condind)
            for i = 1:size(YY, 2)
                Cy{c, i} = squeeze(sum(cat(3, YY{condind{c}, i}), 3))./((nsamples-1)*length(condind{c}));
            end
        end
    case 'bytrial'
        for c = 1:ntrials
            for i = 1:size(YY, 2)
                Cy{c, i} = YY{c, i}./(nsamples-1);
            end
        end
        
end

spm('Pointer', 'Watch');drawnow;

scale = eye(nchan);

if S.scale == 1
    
    sigma = svd(sumYY);
    
    disp('Using spm_pca_order to get scaling factor');
    [M_opt,log_ev] = spm_pca_order(sumYY, N);
    
    scale = scale.*sum(sigma(M_opt:end))./nchan;         
end

for c = 1:size(Cy, 1)
    spm_progress_bar('Init', nvert, ...
        sprintf('Scanning grid points image %d/%d', c, size(Cy, 1))); drawnow;
    if nvert > 100, Ibar = floor(linspace(1, nvert,100));
    else Ibar = 1:nvert; end
    
    pow = nan(1, nvert);
    cpow = zeros(1, numel(Cy(c, :)));
    
    for i = 1:nvert
        if ~isnan(W{i})
            w    = W{i};
            
            for j = 1:numel(Cy(c, :))
                p = real(w*Cy{c, j}*w');
                if isscalar(p) || isequal(S.powermethod, 'trace')
                    cpow(j) = trace(p);
                    
                    if S.scale
                        np = trace(w*scale*w');
                        cpow(j) = cpow(j)./np;
                    end
                elseif isequal(S.powermethod, 'lambda1')
                    [u,s,v] = svd(real(p));
                    cpow(j) = s(1);
                     if S.scale
                        w  = u(:, 1)'*w;
                        np = w*scale*w';
                        cpow(j) = cpow(j)./np;
                    end
                else
                    error('Unsupported option')
                end
            end
            
            if S.logpower
                pow(i) = log(cpow)*S.contrast';
            else
                pow(i) = cpow*S.contrast';
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
            image(c).label = ['uv_pow_'  spm_file(D.fname, 'basename')];
        case 'bycondition'
            image(c).label = ['uv_pow_cond_' S.whatconditions.condlabel{c} '_' spm_file(D.fname, 'basename')];
        case 'bytrial'
            for k = 1:numel(condind)
                if any(c == condind{k})
                    break;
                end
            end
            image(c).label = ['uv_pow_cond_' S.whatconditions.condlabel{k}...
                '_trial_' num2str(alltrials(c)) '_' spm_file(D.fname, 'basename')];
    end
    
end

spm('Pointer', 'Arrow');drawnow;

res = image;