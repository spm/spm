function D = spm_eeg_detect_eyeblinks(S)
% Detects eyeblinks in spm continuous data file
% FORMAT  D = spm_eeg_detect_eyeblinks(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%          .stdthresh  - threshold to reject things that look like
%                         eye-blinks but probably aren't (default: 3)
%          .overwrite  - 1 - replace previous eybelink events (default)
%                        0 - append
% Output:
% D                 - MEEG object with added eyeblink events(also
%                     written on disk)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Laurence Hunt
% $Id: spm_eeg_detect_eyeblinks.m 3858 2010-04-30 20:56:20Z vladimir $

SVNrev = '$Rev: 3858 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Eyeblink detect'); spm('Pointer','Watch');


%-Test for the presence of required Matlab toolbox
%--------------------------------------------------------------------------
if ~license('test','signal_toolbox')
    error('Signal Processing Toolbox is required for eyeblink detection.');
end

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

if ~strcmp(D.type, 'continuous')
    error('Only continuous data is supported at the moment.');
end


% Get the indicex for EOG channel
%--------------------------------------------------------------------------
if ~(isfield(S, 'eogchan') && ~isempty(S.eogchan))
   eogchan = setdiff(D.eogchannels, D.badchannels);
   if length(eogchan)~=1
       [selection, ok]= listdlg('ListString', D.chanlabels(eogchan), 'SelectionMode', 'single' ,'Name', 'Select EOG channel' , 'ListSize', [400 300]);
       if ~ok
           return;
       end
       if ~isempty(eogchan)
           eogchan = eogchan(selection);
       else
           eogchan = selection;
       end
       S.eogchan = D.chanlabels(eogchan);
   end
elseif ~isnumeric(S.eogchan)
    eogchan = D.indchannel(S.eogchan);
else
    eogchan = S.eogchan;
end

try 
    stdthresh = S.stdthresh;
catch
    stdthresh = 3;
end

if ~isfield(S, 'overwrite')
    S.overwrite = spm_input('Overwrite previous?','+1','yes|no',[1 0], 1);
end

%% get EOG data
if length(eogchan)~=1
    error('More than one EOG channel - not currently supported')
end

eog_data = D(eogchan,:,:);

%% filter data at 1-15Hz (eyeblink duration typically 100-300ms) and demean
eog_filt = detrend(ft_preproc_bandpassfilter(eog_data, D.fsample, [1 15], 1001, 'fir'), 'constant');

%% find eye-movements

sd_eeg = std(eog_filt);
em_thresh = stdthresh*sd_eeg;

%% find 'spikes' (putative eyeblinks):

eblength = round(D.fsample/5); %length of eyeblink(200 ms) in samples;
spikes = [];
for i = eblength:length(eog_filt)-eblength;
    if abs(eog_filt(i))>em_thresh && ... %bigger than threshold
       all(abs(eog_filt(i))>=abs(eog_filt(i-eblength+1:i+eblength))); %biggest in surrounding 400ms
        spikes = [spikes i];
    end
end

if isempty(spikes)
    error('No eye-blinks detected by algorithm. Try a lower threshold.')
end

spikemat = zeros(eblength*2,length(spikes));
for i = 1:length(spikes)
    spikemat(:,i) = eog_filt(spikes(i)-eblength+1:spikes(i)+eblength);
end

%reject spikes whose peak is not within 0.5 s.d. of the mean (gets rid of most artefacts
%    etc. not removed by filtering):
mn_spike = mean(spikemat(eblength,:));
sd_spike = 0.5*std(spikemat(eblength,:));
spikes(spikemat(eblength,:)>mn_spike+sd_spike | ...
       spikemat(eblength,:)<mn_spike-sd_spike) = [];
spikemat(:,find(spikemat(eblength,:)>mn_spike+sd_spike | ...
       spikemat(eblength,:)<mn_spike-sd_spike)) = [];

disp(['Number of putative eyeblinks detected: ' num2str(length(spikes))]);
          
% plot
%----------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
colormap(gray)
figure(Fgraph)
clf
subplot(2, 1 , 1)
plot(spikes,ones(length(spikes),1)*5*sd_eeg,'r.');
hold on;
plot(eog_filt);

subplot(2, 1 , 2)
hold on;
plot(spikemat);plot(mean(spikemat,2),'Color','k','LineWidth',4);


% Update the event structure
%----------------------------------------------------------------------
if ~isempty(spikes)  
    for n = 1:D.ntrials
        cspikes   = spikes(spikes>(D.nsamples*(n-1)) & spikes<(D.nsamples*n));
        ctime  = D.trialonset(n)+(cspikes - D.nsamples*(n-1)-1)/D.fsample;
        ctime  = num2cell(ctime);
        
        ev = events(D, n);
        
        if iscell(ev)
            ev = ev{1};
        end
        
        
        if ~isempty(ev) && S.overwrite
            ind1 = strmatch('artefact', {ev.type}, 'exact');
            if ~isempty(ind1)
                ind2 = strmatch('eyeblink', {ev(ind1).value}, 'exact');
                if ~isempty(ind2)
                    ev(ind1(ind2)) = [];
                end
            end
        end
        
        Nevents = numel(ev);
        for i=1:numel(ctime)
            ev(Nevents+i).type     = 'artefact';
            ev(Nevents+i).value    = 'eyeblink';
            ev(Nevents+i).duration = [];
            ev(Nevents+i).time     = ctime{i};
        end
        
        if ~isempty(ev)
            [tevent, I] = sort([ev.time]);
            ev = ev(I);
            D = events(D, n, ev);
        end
    end    
else
    warning(['No eye blinks events detected in the selected channel']);
end

%-Update history (not necessary; leave as call to spm_eeg_montage?) Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
D = D.history(mfilename, S);

save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','Eyeblink detect: done'); spm('Pointer','Arrow');




    
