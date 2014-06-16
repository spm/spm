function res = spm_eeg_artefact_eyeblink(S)
% Detects eyeblinks in spm continuous data file
% S                     - input structure
% fields of S:
%    S.D                - M/EEG object
%    S.chanind          - vector of indices of channels that this plugin will look at.
%    S.threshold        - threshold parameter (in stdev)
%
%    Additional parameters can be defined specific for each plugin
% Output:
%  res -
%   If no input is provided the plugin returns a cfg branch for itself
%
%   If input is provided the plugin returns a matrix of size D.nchannels x D.ntrials
%   with zeros for clean channel/trials and ones for artefacts.
%______________________________________________________________________________________
% Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging

% Laurence Hunt
% $Id: spm_eeg_artefact_eyeblink.m 6046 2014-06-16 10:58:27Z vladimir $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_artefact
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0
    threshold = cfg_entry;
    threshold.tag = 'threshold';
    threshold.name = 'Threshold';
    threshold.strtype = 'r';
    threshold.val = {4};
    threshold.num = [1 1];
    threshold.help = {'Threshold to reject things that look like eye-blinks but probably aren''t'};
           
    excwin = cfg_entry;
    excwin.tag = 'excwin';
    excwin.name = 'Excision window';
    excwin.strtype = 'r';
    excwin.num = [1 1];
    excwin.val = {0};
    excwin.help = {'Window (in ms) to mark as bad around each eyeblink, 0 to not mark data as bad'};
    
    eyeblink = cfg_branch;
    eyeblink.tag = 'eyeblink';
    eyeblink.name = 'Eyeblinks';
    eyeblink.val = {threshold, excwin};
    
    res = eyeblink;
    
    return
end

SVNrev = '$Rev: 6046 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG eyeblink detection');

if isequal(S.mode, 'reject')
     error('Only mark mode is supported by this plug-in, use event-based rejection to reject.');
end

D = spm_eeg_load(S.D);

chanind  =  S.chanind;

if length(chanind)~=1
    error('More than one channel - not currently supported')
end

eog_data = reshape(squeeze(D(chanind,:,:)), 1, []);

%% filter data at 1-15Hz (eyeblink duration typically 100-300ms) and demean
eog_filt = detrend(abs(hilbert(ft_preproc_bandpassfilter(eog_data, D.fsample, [1 15], 4, 'but',[], 'split'))), 'constant');

%% find eye-movements

sd_eeg=(spm_percentile(eog_filt,85)-spm_percentile(eog_filt,15))/2; %robust estimate of standard deviation, suggested by Mark Woolrich
em_thresh = S.threshold*sd_eeg;

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

%reject spikes whose peak is not within 1 s.d. of the mean (gets rid of most artefacts
%    etc. not removed by filtering):
mn_spike = mean(spikemat(eblength,:));
sd_spike = std(spikemat(eblength,:));
spikes(spikemat(eblength,:)>mn_spike+sd_spike | ...
       spikemat(eblength,:)<mn_spike-sd_spike) = [];
spikemat(:,find(spikemat(eblength,:)>mn_spike+sd_spike | ...
       spikemat(eblength,:)<mn_spike-sd_spike)) = [];

disp(['Number of putative eyeblinks detected: ' num2str(length(spikes))]);
       
num_eb_per_min=(60*D.fsample*length(spikes))/length(eog_data);
disp([num2str(num_eb_per_min) ' eye-blinks per minute '])
if (num_eb_per_min<0.5)
    error(['Only ' num2str(num_eb_per_min) ' eye-blinks per minute detected by algorithm. Try a lower threshold.'])
end
if (num_eb_per_min>60)
    error(['As many as ' num2str(num_eb_per_min) ' eye-blinks per minute detected by algorithm. Try a higher threshold.'])
end

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
        
        
        if ~isempty(ev) && ~S.append
            ind1 = strmatch('artefact_eyeblink', {ev.type}, 'exact');
            if ~isempty(ind1)
                ind2 = strmatch(D.chanlabels(chanind), {ev(ind1).value}, 'exact');
                if ~isempty(ind2)
                    ev(ind1(ind2)) = [];
                end
            end
        end        
        
        Nevents = numel(ev);
        for i=1:numel(ctime)
            if ctime{i} == 0
                continue; %likely to be trial border falsely detected as eyeblink
            end
            ev(Nevents+i).type     = 'artefact_eyeblink';
            ev(Nevents+i).value    = char(D.chanlabels(chanind));
            if S.excwin == 0
                ev(Nevents+i).duration = [];
                ev(Nevents+i).time     = ctime{i};
            else
                ev(Nevents+i).time     = max(D.trialonset(n), ctime{i} - 5e-4*S.excwin);
                ev(Nevents+i).duration = min(1e-3*S.excwin, (D.time(end)-D.time(1))-(ev(Nevents+i).time-D.trialonset(n)))+...
                    min(ctime{i} - 5e-4*S.excwin, 0);               
            end                      
        end
        
        if ~isempty(ev)
            [dum, I] = sort([ev.time]);
            ev = ev(I);
            D = events(D, n, ev);
        end
    end    
else
    warning(['No eye blinks events detected in the selected channel']);
end

res = D;

spm('FigName','M/EEG eyeblink detection: done');