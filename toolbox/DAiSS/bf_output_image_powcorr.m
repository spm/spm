function res = bf_output_image_powcorr(BF, S)
% Computes phase-amplitude coupling
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak using code bits from OSL library
% $Id: bf_output_image_powcorr.m 7703 2019-11-22 12:06:29Z guillaume $

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
    woi.name = 'Time window of interest';
    woi.strtype = 'r';
    woi.num = [1 2];
    woi.val = {[-Inf Inf]};
    woi.help = {'Time windows (in ms)'};
    
    freqref         = cfg_entry;
    freqref.tag     = 'freqref';
    freqref.name    = 'Reference frequencies';
    freqref.strtype = 'r';
    freqref.num     = [1 Inf];
    freqref.val     = {20};
    freqref.help    = {'Frequencies in the reference channel'};
    
    resref         = cfg_entry;
    resref.tag     = 'resref';
    resref.name    = 'Reference resolutions';
    resref.strtype = 'r';
    resref.num     = [1 Inf];
    resref.val     = {5};
    resref.help    = {'Frequency resolution for reference frequencies. Single value or vector'};
    
    freq         = cfg_entry;
    freq.tag     = 'freq';
    freq.name    = 'Data frequencies';
    freq.strtype = 'r';
    freq.num     = [1 Inf];
    freq.val     = {20};
    freq.help    = {'First set of frequencies'};
    
    res         = cfg_entry;
    res.tag     = 'res';
    res.name    = 'Data resolutions';
    res.strtype = 'r';
    res.num     = [1 Inf];
    res.val     = {5};
    res.help    = {'Frequency resolution for data frequencies. Single value or vector'};
    
    refchan = cfg_entry;
    refchan.tag = 'refchan';
    refchan.name = 'Reference channel';
    refchan.strtype = 's';
    refchan.num = [1 Inf];
    refchan.help = {'Reference channel name.'};
    
    movavg         = cfg_entry;
    movavg.tag     = 'movavg';
    movavg.name    = 'Moving average window';
    movavg.strtype = 'r';
    movavg.num     = [1 1];
    movavg.val     = {100};
    movavg.help    = {'Time window for moving average of power envelope (ms).',...
        'Specify 0 to not average'};
    
    movavg         = cfg_entry;
    movavg.tag     = 'movavg';
    movavg.name    = 'Moving average window';
    movavg.strtype = 'r';
    movavg.num     = [1 1];
    movavg.val     = {100};
    movavg.help    = {'Time window for moving average of power envelope (ms).',...
        'Specify 0 to not average'};
    
    subsample = cfg_entry;
    subsample.tag = 'subsample';
    subsample.name = 'Subsample';
    subsample.strtype = 'n';
    subsample.num = [1 1];
    subsample.val = {1};
    subsample.help = {'Set to N to subsample the power to every Nth sample'};
    
    shuffle         = cfg_entry;
    shuffle.tag     = 'shuffle';
    shuffle.name    = 'Shuffle';
    shuffle.strtype = 'w';
    shuffle.num     = [1 1];
    shuffle.help    = {'Shuffle the reference channel to produce the null case.',...
        'Specify the number of shufflings'};
    shuffle.val = {0};
    
    
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
    
    image_powcorr      = cfg_branch;
    image_powcorr.tag  = 'image_powcorr';
    image_powcorr.name = 'Power correlations image';
    image_powcorr.val  = {whatconditions, sametrials, shuffle, woi, refchan, freqref, ....
        resref, freq, res, movavg, subsample, modality};
    
    res = image_powcorr;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

S.woi = 1e-3*S.woi; % ms -> s

samples =  D.indsample(S.woi(1)):D.indsample(S.woi(2));
nsamples = length(samples);
times = D.time(samples);

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

alltrials = spm_vec(trials);
ntrials   = length(alltrials);

nref = length(S.freqref);
nfreq    = length(S.freq);

W = BF.inverse.(S.modality).W;
nvert = numel(W);

Y = U'*reshape(D(channels, samples, alltrials), nchan, []);
Y = reshape(Y, size(Y, 1), nsamples, ntrials);

Yr = squeeze(D(D.indchannel(S.refchan), samples, alltrials));
if size(Yr, 1) == 1
    Yr = Yr';
end

spectrum = ft_specest_hilbert(Yr', times,...
    'freqoi', S.freqref, 'width', S.resref, 'filttype', 'but', 'filtorder', 2,...
    'filtdir', 'twopass', 'verbose', 0);

spectrum = reshape(permute(spectrum, [3 1 2]), nsamples, ntrials*nref);
spectrum = abs(spectrum);

if S.movavg
    avwin =  spm_hanning(1e-3*(S.movavg)*D.fsample);
    spectrum = conv2(avwin, 1, spectrum, 'same');
end

spectrum = spectrum(1:S.subsample:end, :);

spectrum = detrend(spectrum);
spectrum = spectrum./repmat(std(spectrum), size(spectrum, 1), 1);

refsig   = reshape(spectrum, [], ntrials, nref);

powc   = nan(length(S.freqref), length(S.freq), nvert);

if S.shuffle
    spowc = repmat(powc, [1 1 1 S.shuffle]);
    sind = zeros(S.shuffle, ntrials);
    for s = 1:S.shuffle
        sind(s, :) = randperm(ntrials);
    end
end

for f = 1:length(S.freq)
    
    spm_progress_bar('Init', ntrials, ...
        sprintf('Computing data spectra')); drawnow;
    if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
    else Ibar = 1:ntrials; end
    
    
    Yh = 0*Y;
    for i = 1:ntrials
        Yh(: , : ,i) = spm_squeeze(ft_specest_hilbert(squeeze(Y(:,:, i)), times,...
            'freqoi', S.freq(f), 'width', S.res(f), 'filttype', 'but', ...
            'filtorder', 2,  'filtdir', 'twopass', 'verbose', 0), 2);
        
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end
    end
    
    Yh = reshape(Yh, nchan, []);
    
    spm_progress_bar('Clear');
    
    spm_progress_bar('Init', nvert, ...
        sprintf('Scanning grid points image')); drawnow;
    if nvert > 100, Ibar = floor(linspace(1, nvert,100));
    else Ibar = 1:nvert; end
    
    for i = 1:nvert
        if ~isnan(W{i})
            w    = W{i};
            
            sYh  = w*Yh;
            
            sYh  = reshape(abs(sYh), nsamples, ntrials);
            
            if S.movavg
                sYh = conv2(avwin, 1, sYh, 'same');
            end
            
            sYh = sYh(1:S.subsample:end, :);
            
            sYh = detrend(sYh);
            
            sYh = sYh./repmat(std(sYh), size(sYh, 1), 1);
            
            x = sYh(:);
            
            for j = 1:nref
                
                rYh = spm_squeeze(refsig(:,:, j), 3);
                
                for shuffle = 0:S.shuffle
                    if shuffle
                        rYh = rYh(:, sind(shuffle, :));
                    end
                    
                    y = rYh(:);
                    
                    pinvx = pinv(x);
                    pe  = pinvx*y;
                    r   = y-x*pe;
                    vr  = diag(r'*r/(size(y,1)-size(x,2)));
                    vrp = pinv(x'*x)*vr;
                    cs  = pe/sqrt(vrp);
                    
                    if shuffle
                        spowc(j, f, i, shuffle) = cs;
                    else
                        powc(j, f, i) = cs;
                    end
                end
            end
        end
        
        
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end
    end
end

spm_progress_bar('Clear');



if max(nfreq, nref)>1
    image(1).val     = squeeze(sum(sum(powc, 2), 1));
    image(1).label   = ['sumpowc_'  spm_file(D.fname, 'basename')];
    c = 2;
else
    c = 1;
end

for f = 1:nref
    for g = 1:nfreq
        image(c).val     = squeeze(powc(f, g, :));
        image(c).label   = ['powc_ref_' num2str(S.freqref(f)) 'Hz_meg_'...
            num2str(S.freq(g)) 'Hz_' spm_file(D.fname, 'basename')];
        c = c+1;
    end
end

for shuffle = 1:S.shuffle
    if max(nfreq, nref)>1
        image(c).val     = squeeze(sum(sum(spowc(:,:,:, shuffle), 2), 1));
        image(c).label   = ['shuffled' num2str(shuffle) '_sumpowc_'  spm_file(D.fname, 'basename')];
        c = c+1;
    end
    
    for f = 1:nref
        for g = 1:nfreq
            image(c).val     = squeeze(spowc(f, g, :, shuffle));
            image(c).label   = ['shuffled' num2str(shuffle) '_powc_ref_' num2str(S.freqref(f)) 'Hz_meg_'...
                num2str(S.freq(g)) 'Hz_' spm_file(D.fname, 'basename')];
            c = c+1;
        end
    end
end

res = image;