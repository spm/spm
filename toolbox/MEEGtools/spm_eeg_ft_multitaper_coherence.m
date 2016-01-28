function Dcoh = spm_eeg_ft_multitaper_coherence(S)
% Function for computing time-frequency decomposition using multitaper
%
% WARNING: This function uses some quite specific settings and is not generic. It is
% just an example of how Fieldtrip spectral analysis can be combined with
% SPM
%
% FORMAT D = spm_eeg_ft_multitaper_coherence(S)
%
% S           - input structure (optional)
% (optional) fields of S:
%   S.D        - filename, or M/EEG object
%   S.chancomb - Nx2 cell array with channel pairs
%   S.pretrig  - time to start TF analysis in PST (ms)
%   S.posttrig - time to end TF analysis in PST (ms)
%   S.timewin  - time window (resolution) in ms
%   S.timestep - time step in ms
%   S.freqwin  - frequency window (Hz)
%   S.freqres  - frequency resolution
%   S.robust      - (optional) - use robust averaging for computing
%                                coherence
%                 .savew  - save the weights in an additional dataset
%                 .bycondition - compute the weights by condition (1,
%                                default) or from all trials (0)
%                 .ks     - offset of the weighting function (default: 3)
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Institute of Neurology, UCL

% Vladimir Litvak
% $Id: spm_eeg_ft_multitaper_coherence.m 6699 2016-01-28 09:55:26Z vladimir $

%%
SVNrev = '$Rev: 6699 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Fieldtrip multitaper coherence');

%%
%-Test for the presence of required Matlab toolbox
%--------------------------------------------------------------------------
if ~license('test','signal_toolbox')
    error('This function requires the Signal Processing Toolbox.');
end

%% ============ Load SPM EEG file and verify consistency
if nargin == 0
    S = [];
end

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
    S.D = D;
end

if ischar(D)
    try
        D = spm_eeg_load(D);
    catch
        error(sprintf('Trouble reading file %s', D));
    end
end

%-Configure the spectral analysis
%--------------------------------------------------------------------------
if ~isfield(S, 'chancomb')
    switch spm_input('What channel combinations', 1,'all|meeg|select', char('all', 'meeg', 'select'), 'all')
        case 'all'
            S.chancomb  = ft_channelcombination('all', D.chanlabels);
        case 'meeg'
            S.chancomb = ft_channelcombination('all', D.chanlabels(D.indchantype('MEEG', 'GOOD')));
        case 'select'
            S.chancomb = {};
            while 1
                [selection, ok]= listdlg('ListString', D.chanlabels, 'SelectionMode', 'multiple' ,'Name', 'Select two channels' , 'ListSize', [400 300]);
                if ~ok
                    break;
                elseif length(selection)==2
                    S.chancomb = [S.chancomb; D.chanlabels(selection)];
                end
            end
    end
end

S.chancomb = S.chancomb(all(ismember(S.chancomb, D.chanlabels)'), :);

if ~isfield(S, 'pretrig')
    S.pretrig = spm_input('Start in PST [ms]', '+1', 'r', '', 1);
end

if ~isfield(S, 'posttrig')
    S.posttrig = spm_input('End in PST [ms]', '+1', 'r', '', 1);
end

if ~isfield(S, 'timewin')
    S.timewin = spm_input('Time window (ms)', '+1', 'r', '400', 1);
end

if ~isfield(S, 'timestep')
    S.timestep = spm_input('Time step (ms)', '+1', 'r', '50', 1);
end

if ~isfield(S, 'freqwin')
    S.freqwin = spm_input('Frequency window (Hz)', '+1', 'r', '0 90', 2);
end

if ~isfield(S, 'robust') && spm_input('Use robust averaging?','+1','yes|no', [1 0], 0)
    S.robust = [];
else
    S.robust = 'no';
end
        
if  ~isequal(S.robust, 'no')
    robust = 1;
    try
        ks =  S.robust.ks;
    catch
        ks = spm_input('Offset weighting function by', '+1', 'r', '3', 1);
        S.robust.ks = ks;
    end

    if ~isfield(S.robust, 'savew')
        savew =  spm_input('Save weights?','+1','yes|no', [1 0], 1);
        S.robust.savew = savew;
    else
        savew = S.robust.savew;
    end

    if ~isfield(S.robust, 'bycondition')
        bycondition =  spm_input('Compute weights by condition?','+1','yes|no', [1 0], 1);
        S.robust.bycondition = bycondition;
    else
        bycondition = S.robust.bycondition;
    end
else
    robust = 0;
end

spm('Pointer','Watch');

data = D.ftraw(D.indchantype('ALL', 'GOOD'), ':', ':');


prestim = 1e-3*S.pretrig;
poststim = 1e-3*S.posttrig;
%%
timewin = 1e-3*S.timewin;
step = 1e-3*S.timestep;

%-Run the Fieldtrip code
%--------------------------------------------------------------------------

cfg = [];
cfg.output ='powandcsd';
cfg.taper = 'dpss';
cfg.channel = unique(S.chancomb(:));
cfg.channelcmb = S.chancomb;
cfg.method          = 'mtmconvol';
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';

% This sets the centers of frequency bins at the optimal locations based on
% the time window.
cfg.foi             = (1/timewin):(1/timewin):S.freqwin(2); % Frequency axis
cfg.foi             = cfg.foi(cfg.foi>=S.freqwin(1));
numfoi              = length(cfg.foi);

% This means that the time resolution is the same for all frequencies
cfg.t_ftimwin       = zeros(1,numfoi);
cfg.t_ftimwin(:)    = timewin; % Time resolution

% This part is about frequency resolution
cfg.tapsmofrq       = zeros(1,numfoi); % Initialize to zero
cfg.tapsmofrq(:)    = 1/timewin; % Set initial resolution to 1/timewin (i.e. 2.5 Hz) for all frequencis
% Here it sets the resolution for frequencies above 10*(1/timewin) (25 Hz)
% to 0.1 times the frequency. This means that at up to 25 Hz the resolution
% is fixed and then it starts slowly increasing up to 10Hz in each
% direction for 100 Hz. If you comment out this line, you'll have fixed
% time resolution for all frequencies.
cfg.tapsmofrq(cfg.foi>10*(1/timewin))    = 0.1*cfg.foi(cfg.foi>10*(1/timewin));
cfg.tapsmofrq(cfg.foi>50)                = 5;

if isfield(S, 'freqres')
    cfg.tapsmofrq(:)                         = S.freqres;
end
% figure; plot(cfg.foi, cfg.tapsmofrq);xlabel('frequency (Hz)');ylabel('frequency resolution (Hz)')
% This is the time axis. The window with the width defined above (400 msec)
% is moved every time by 'step' (100 ms). The earliest you can start is
% half the time window from the start of the data. Otherwise your time
% window will overlap a segment with no data and you will get NaNs in the
% output. The same idea at the end.
cfg.toi=(prestim+(timewin/2)):step:(poststim-(timewin/2)-1/D.fsample); % Time axis

freq = ft_freqanalysis(cfg, data);

if ~isfield(freq, 'time')
    freq.time = cfg.toi;
    freq.dimord = [freq.dimord '_time'];
end

S.chancomb = ft_channelcombination(S.chancomb, freq.label);

np = size(S.chancomb, 1);

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
Dcoh = clone(D, ['COH' fname(D)], [np length(freq.freq) length(freq.time) D.nconditions]);
Dcoh = frequencies(Dcoh, ':', freq.freq);
Dcoh = fsample(Dcoh, 1./mean(diff(freq.time)));
Dcoh = timeonset(Dcoh, freq.time(1));
Dcoh = check(Dcoh);

if robust && savew
    Dw = clone(Dcoh, ['WCOH' fnamedat(D)], [np length(freq.freq) 3*length(freq.time) D.ntrials]);
    Dw = frequencies(Dw, freq.freq);
    Dw = fsample(Dw, 1./mean(diff(freq.time)));
    Dw = timeonset(Dw, freq.time(1));
    Dw = check(Dw);
end
%%
%--------------------------------------------------------------------------
cl = D.condlist;
nc = length(cl);

spm_progress_bar('Init', np, 'Pairs completed');

ni = zeros(1,D.nconditions);
for i = 1:D.nconditions
    w = indtrial(D, deblank(cl{i}), 'GOOD')';
    ni(i) = length(w);
    if ni(i) == 0
        warning('%s: No trials for trial type %d', D.fname, cl{i});
    end
end

goodtrials  =  indtrial(D, cl, 'GOOD');

for j = 1:np
    powind   = find(ismember(freq.label, S.chancomb(j, :)));
    crossind = find(all(ismember(freq.labelcmb, S.chancomb(j, :))'));

    if robust && ~bycondition
        w = goodtrials;
        [Y, W1] = spm_robust_average(permute(cat(4, spm_squeeze(abs(freq.crsspctrm(w, crossind , :, :)), 2),...
            spm_squeeze(freq.powspctrm(w, powind(1), :, :), 2), spm_squeeze(freq.powspctrm(w, powind(2), :, :), 2)), ...
            [2 3 1 4]), 3, ks);

        if savew
            Dw(j, :, :, goodtrials) = ...
                shiftdim(cat(2, W1(:, :, :, 1), W1(:, :, :, 2), W1(:, :, :, 3)), -1);
        end
        W = zeros([length(freq.freq) length(freq.time) D.ntrials 3]);
        W(:, :, goodtrials, :) = W1;
    end

    for i = 1:nc
        w = D.indtrial(cl{i}, 'GOOD');

        if length(w)<2
            continue;
        end

        if ~robust
            Dcoh(j, :, :, i) = spm_squeeze(abs(mean(freq.crsspctrm(w, crossind, :, :))./...
                sqrt(mean(freq.powspctrm(w, powind(1), :, :)).*mean(freq.powspctrm(w, powind(2), :, :)))), 2);
        else
            if bycondition
                [Y, W] = spm_robust_average(permute(cat(4, spm_squeeze(abs(freq.crsspctrm(w, crossind , :, :)), 2),...
                    spm_squeeze(freq.powspctrm(w, powind(1), :, :), 2), spm_squeeze(freq.powspctrm(w, powind(2), :, :), 2)), ...
                    [2 3 1 4]), 3, ks);
                
                cross = spm_squeeze(permute(freq.crsspctrm(w, crossind, :, :), [2 3 4 1]), 1);
                WW =  spm_squeeze(W(:, :, :, 1), 4);
                cross(WW == 0) = 0; % This is to get rid of NaNs because NaN*0 == NaN
                cross = abs(sum(cross.*WW, 3))./spm_squeeze(sum(WW, 3), 3);
                
                pow1 = spm_squeeze(permute(freq.powspctrm(w, powind(1), :, :), [2 3 4 1]), 1);
                WW = spm_squeeze(W(:, :, :, 2), 4);
                pow1(WW == 0) = 0;
                pow1 = sum(pow1.*WW, 3)./spm_squeeze(sum(WW, 3), 3);
                
                pow2 = spm_squeeze(permute(freq.powspctrm(w, powind(2), :, :), [2 3 4 1]), 1);
                WW = spm_squeeze(W(:, :, :, 3), 4);
                pow2(WW == 0) = 0;
                pow2 = sum(pow2.*WW, 3)./spm_squeeze(sum(WW, 3), 3);
                               
                Dcoh(j, :, :, i) = shiftdim(cross./sqrt(pow1.*pow2), -1);

                if savew
                    Dw(j, :, :, w) = ...
                        shiftdim(cat(2, W(:, :, :, 1), W(:, :, :, 2), W(:, :, :, 3)), -1);
                end
            else
                
                cross = spm_squeeze(permute(freq.crsspctrm(w, crossind, :, :), [2 3 4 1]), 1);
                WW =  spm_squeeze(W(:, :, w, 1), 4);
                cross(WW == 0) = 0;
                cross = abs(sum(cross.*WW, 3))./spm_squeeze(sum(WW, 3), 3);
                
                pow1 = spm_squeeze(permute(freq.powspctrm(w, powind(1), :, :), [2 3 4 1]), 1);
                WW = spm_squeeze(W(:, :, w, 2), 4);
                pow1(WW == 0) = 0;
                pow1 = sum(pow1.*WW, 3)./spm_squeeze(sum(WW, 3), 3);
                
                pow2 = spm_squeeze(permute(freq.powspctrm(w, powind(2), :, :), [2 3 4 1]), 1);
                WW = spm_squeeze(W(:, :, w, 3), 4);
                pow2(WW == 0) = 0;
                pow2 = sum(pow2.*WW, 3)./spm_squeeze(sum(WW, 3), 3);               
                
                Dcoh(j, :, :, i) = shiftdim(cross./sqrt(pow1.*pow2), -1);
            end
        end
    end
    
    Dcoh = chanlabels(Dcoh, j, [S.chancomb{j, 1} ' vs. ' S.chancomb{j, 2}]);
    
    if isequal(D.chantype(strmatch(S.chancomb{j, 1}, D.chanlabels, 'exact')), ...
            D.chantype(strmatch(S.chancomb{j, 2}, D.chanlabels, 'exact')))
        Dcoh = chantype(Dcoh, j, D.chantype(strmatch(S.chancomb{j, 1}, D.chanlabels, 'exact')));
    end      
    
    
    spm_progress_bar('Set', j);
end
%%

%-Copy some additional information from the original file
%--------------------------------------------------------------------------
Dcoh  = conditions (Dcoh, ':', cl);
Dcoh  = repl(Dcoh, ':', ni);

Dcoh = history(Dcoh, history(D));

%-Update history
%--------------------------------------------------------------------------
Dcoh = history(Dcoh, mfilename, S);

%-Save
%--------------------------------------------------------------------------
save(Dcoh);

if robust && savew
    Dw = chanlabels(Dw, [], chanlabels(Dcoh));
    Dw = chantype(Dw, [], chantype(Dcoh));
    save(Dw);
end
    
%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','M/EEG merge: done'); spm('Pointer','Arrow');


