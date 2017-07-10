function res = spm_eeg_artefact_threshchan(S)
% Plugin for spm_eeg_artefact doing artefact detection by chanel thresholding
% S            - input structure
% fields of S:
%    S.D       - M/EEG object
%    S.chanind - vector of indices of channels that this plugin will look at
%                         
%    Additional parameters can be defined specific for each plugin.
%
% Output:
% res - 
%    If no input is provided the plugin returns a cfg branch for itself.
%
%    If input is provided the plugin returns a matrix of size D.nchannels x D.ntrials  
%    with zeros for clean channel/trials and ones for artefacts.
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_artefact_threshchan.m 7132 2017-07-10 16:22:58Z guillaume $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_artefact
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0
    threshold         = cfg_entry;
    threshold.tag     = 'threshold';
    threshold.name    = 'Threshold';
    threshold.strtype = 'r';
    threshold.num     = [1 1];
    threshold.help    = {'Threshold value to apply to all channels.'};

    excwin         = cfg_entry;
    excwin.tag     = 'excwin';
    excwin.name    = 'Excision window';
    excwin.strtype = 'r';
    excwin.num     = [1 1];
    excwin.val     = {1000};
    excwin.help    = {'Window (in ms) to mark as bad around each jump (for mark mode only), 0 - do not mark data as bad.'};
    
    
    threshchan      = cfg_branch;
    threshchan.tag  = 'threshchan';
    threshchan.name = 'Threshold channels';
    threshchan.val  = {threshold, excwin};
    threshchan.help = {''};
    
    res = threshchan;
    
    return
end

SVNrev = '$Rev: 7132 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG threshold channels');

D = spm_eeg_load(S.D);

chanind = S.chanind;
threshold = S.threshold;
res = zeros(D.nchannels, D.ntrials);

if isequal(S.mode, 'reject')
    
    %-Artefact detection
    %----------------------------------------------------------------------
    
    spm_progress_bar('Init', D.ntrials, 'Trials checked');
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
    else Ibar = [1:D.ntrials]; end
    
    for i = 1:D.ntrials
        res(chanind, i) = squeeze(max(abs(D(chanind, :, i)), [], 2))>threshold;
        if any(Ibar == i), spm_progress_bar('Set', i); end
    end
    
    spm_progress_bar('Clear');
    
elseif isequal(S.mode, 'mark')
    multitrial = D.ntrials>1;
    
    if multitrial
        spm_progress_bar('Init', D.ntrials, 'Trials checked');
        if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
        else Ibar = [1:D.ntrials]; end
    end
    
    for i = 1:D.ntrials
        
        bad = abs(squeeze(D(chanind, :, i)))>threshold;
        
        if ~any(bad(:))
            if multitrial && any(Ibar == i), spm_progress_bar('Set', i); end
            continue;
        end
        
        if S.excwin>0
            excwin = ones(1, round(5e-4*S.excwin*D.fsample));
        end
        
        if ~multitrial
            spm_progress_bar('Init', D.nchannels, 'Channels checked');
            if D.nchannels > 100, Ibar = floor(linspace(1, D.nchannels,100));
            else Ibar = [1:D.nchannels]; end
        end
        
        res = [];
        for j = 1:length(chanind)
            bad(j, :) = ~~conv(double(bad(j, :)), excwin, 'same');
         
            if (bad(j, :)/size(bad,2))>S.badchanthresh
                onsets = 1;
                offsets = size(bad,2)+1;
            else
                onsets  = find(bad(j, :));
                offsets = find(~bad(j, :));
                onsets(find(diff(onsets)<2)+1) = [];
                
                if bad(j, end)
                    offsets(end+1) = length(bad)+1;
                end
            end
            
            for k = 1:length(onsets)
                res(end+1).type   = 'artefact_threshold';
                res(end).value    = char(D.chanlabels(chanind(j)));
                res(end).time     = D.time(onsets(k)+1) - D.time(1) + D.trialonset(i);
                res(end).duration = (min(offsets(offsets>onsets(k)))-onsets(k))./D.fsample;
            end
            if ~multitrial && any(Ibar == j), spm_progress_bar('Set', j); end
        end
        
        if ~multitrial, spm_progress_bar('Clear'); end
        
        if ~isempty(res)
            ev = D.events(i);
            if iscell(ev)
                ev = ev{1};
            end
            
            if ~S.append
                ev(strmatch('artefact_threshold', {ev.type})) = [];
            end
            
            D = events(D, i, spm_cat_struct(ev, res));
        end
        
        
        if multitrial && any(Ibar == i), spm_progress_bar('Set', i); end
    end
    
    if multitrial, spm_progress_bar('Clear'); end
    
    res = D;
end
    
spm('FigName', 'M/EEG threshold channels: done');
