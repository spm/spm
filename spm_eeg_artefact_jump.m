function res = spm_eeg_artefact_jump(S)
% Plugin for spm_eeg_artefact doing jump detection.
% S                     - input structure
% fields of S:
%    S.D                - M/EEG object
%    S.chanind          - vector of indices of channels that this plugin will look at.
%                         
%    Additional parameters can be defined specific for each plugin
% Output:
%  res - 
%   If no input is provided the plugin returns a cfg branch for itself
%
%   If input is provided the plugin returns a matrix of size D.nchannels x D.ntrials 
%   with zeros for clean channel/trials and ones for artefacts.
%______________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_artefact_jump.m 5592 2013-07-24 16:25:55Z vladimir $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_artefact
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0
    threshold = cfg_entry;
    threshold.tag = 'threshold';
    threshold.name = 'Threshold';
    threshold.strtype = 'r';
    threshold.num = [1 1];
    threshold.help = {'Threshold value to apply to all channels'};

    excwin = cfg_entry;
    excwin.tag = 'excwin';
    excwin.name = 'Excision window';
    excwin.strtype = 'r';
    excwin.num = [1 1];
    excwin.val = {1000}; 
    excwin.help = {'Window (in ms) to mark as bad around each jump (for mark mode only), 0 - do not mark data as bad'};
    
    jump = cfg_branch;
    jump.tag = 'jump';
    jump.name = 'Difference between adjacent samples';
    jump.val = {threshold, excwin};
    
    res = jump;
    
    return
end

SVNrev = '$Rev: 5592 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG jump detection');

D = spm_eeg_load(S.D);

chanind  = S.chanind;
threshold = S.threshold;
res = zeros(D.nchannels, D.ntrials);


if isequal(S.mode, 'reject')
    res = zeros(D.nchannels, D.ntrials);
    
    %-Artefact detection
    %--------------------------------------------------------------------------
    
    spm_progress_bar('Init', D.ntrials, 'Trials checked');
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
    else Ibar = [1:D.ntrials]; end
    
    for i = 1:D.ntrials
        res(chanind, i) = max(abs(diff(squeeze(D(chanind, :, i)), [], 2)), [], 2)>threshold;
        if ismember(i, Ibar), spm_progress_bar('Set', i); end
    end
    
    spm_progress_bar('Clear');
    
elseif isequal(S.mode, 'mark')   
    
    if isequal(D.type, 'continuous')
        spm_progress_bar('Init', length(chanind), 'Channels checked');
        if length(chanind) > 100, Ibar = floor(linspace(1, length(chanind),100));
        else Ibar = [1:length(chanind)]; end
    else
        spm_progress_bar('Init', D.ntrials, 'Trials checked');
        if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
        else Ibar = [1:D.ntrials]; end
    end
    
    for i = 1:D.ntrials
        res = [];
        for j = 1:length(chanind)
            onsets  = find(abs(diff(squeeze(D(chanind(j), :, i)), [], 2))>threshold);
            if isempty(onsets)
                if isequal(D.type, 'continuous')
                    if ismember(j, Ibar), spm_progress_bar('Set', j); end
                end
                continue;
            end
            
            excwin = round(5e-4*S.excwin*D.fsample);
            wind   = -excwin:excwin;
            excind = repmat(wind, length(onsets), 1) + repmat(onsets(:), 1, length(wind));
            excind = unique(excind);
            excind(excind<1) = [];
            excind(excind>D.nsamples) = [];
            if   (length(excind)/D.nsamples)<S.badchanthresh
                res(end+1).type   = 'artefact_jump';
                res(end).value    = char(D.chanlabels(chanind(j)));
                res(end).time     = D.trialonset(i);
                res(end).duration = D.time(end) - D.time(1);               
            else                              
                k = 1;
                m = 1;
                boundary = onsets(1)+excwin;
                while 1
                    if m <= length(onsets)                        
                        if (onsets(m) - excwin) <= boundary
                            boundary = onsets(m)+excwin;
                            m = m+1;
                        else                                                        
                            res(end+1).type   = 'artefact_jump';
                            res(end).value    = char(D.chanlabels(chanind(j)));
                            res(end).time     = max(0, D.time(onsets(k)) - 5e-4*S.excwin - D.time(1)) + D.trialonset(i);
                            res(end).duration = min(boundary/D.fsample, D.time(end))-res(end).time+D.trialonset(i);
                            
                            k = m+1;
                            m = k;
                        end
                    else
                        break;
                    end
                end
            end
            if isequal(D.type, 'continuous')
                if ismember(j, Ibar), spm_progress_bar('Set', j); end
            end
        end
        if ~isempty(res)
            ev = D.events(i);
            if iscell(ev)
                ev = ev{1};
            end
            
            if ~S.append
                ev(strmatch('artefact_jump', {ev.type})) = [];
            end
            
            D = events(D, i, spm_cat_struct(ev, res));
        end
        
        if ~isequal(D.type, 'continuous')
            if ismember(i, Ibar), spm_progress_bar('Set', i); end
        end
    end
    
    spm_progress_bar('Clear');
    
    res = D;
end

spm('FigName','M/EEG jump detection: done');