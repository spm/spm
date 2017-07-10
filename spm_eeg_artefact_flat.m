function res = spm_eeg_artefact_flat(S)
% Plugin for spm_eeg_artefact doing flat channel detection
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
% $Id: spm_eeg_artefact_flat.m 7132 2017-07-10 16:22:58Z guillaume $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_artefact
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0
    threshold         = cfg_entry;
    threshold.tag     = 'threshold';
    threshold.name    = 'Threshold';
    threshold.strtype = 'r';
    threshold.val     = {0};
    threshold.num     = [1 1];
    threshold.help    = {'Threshold for difference between adjacent samples.'};
    
    
    seqlength         = cfg_entry;
    seqlength.tag     = 'seqlength';
    seqlength.name    = 'Flat segment length';
    seqlength.strtype = 'r';
    seqlength.num     = [1 1];
    seqlength.val     = {4};
    seqlength.help    = {'Minimal number of adjacent samples with the same value to reject.'};
    
    flat      = cfg_branch;
    flat.tag  = 'flat';
    flat.name = 'Flat segments';
    flat.val  = {threshold, seqlength};
    flat.help = {''};
    
    res = flat;
    
    return;
end

SVNrev = '$Rev: 7132 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG flat data detection');

D = spm_eeg_load(S.D);

chanind   =  S.chanind;
threshold = S.threshold;
seqlength = S.seqlength;

if isequal(S.mode, 'reject')
    res = zeros(D.nchannels, D.ntrials);
    
    %-Artefact detection
    %----------------------------------------------------------------------
    
    spm_progress_bar('Init', D.ntrials, 'Trials checked');
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
    else Ibar = [1:D.ntrials]; end
    
    for i = 1:D.ntrials
        for j = 1:length(chanind)
            tmp  = find(abs(diff(squeeze(D(chanind(j), :, i)), [], 2))>threshold);
            if max(diff([0 tmp D.nsamples]))>=seqlength
                res(chanind(j), i) = 1;
            end
        end
        if any(Ibar == i), spm_progress_bar('Set', i); end
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
            dat  = abs(diff(squeeze(D(chanind(j), :, i)), [], 2))>threshold;
            if  sum(dat)/length(dat)<(1-S.badchanthresh)
                res(end+1).type   = 'artefact_flat';
                res(end).value    = char(D.chanlabels(chanind(j)));
                res(end).time     = D.trialonset(i);
                res(end).duration = D.time(end) - D.time(1) + 1;
            else
                tmp  = find(dat);
                diffs = diff([0 tmp D.nsamples]);
                onsets = find(diffs>=seqlength);
                
                onsetsamples = [];
                if any(onsets == 1);
                    onsetsamples = 1;
                    onsets(1)    = [];
                    onsetsamples = [onsetsamples tmp(onsets-1)+1];
                    onsets       = [1 onsets];
                else
                   onsetsamples = [onsetsamples tmp(onsets-1)+1];
                end                   
                                                   
                k = 1;
                m = 1;
                while k<=length(onsets)
                    if m <= length(onsets)
                        ind1 = onsetsamples(k);
                        ind2 = onsetsamples(m) + diffs(onsets(m))-2;
                        if ind2 > length(dat)
                            ind2 = length(dat);
                        end
                        if (sum(dat(ind1:ind2))/(ind2-ind1+1))<0.5 
                            m = m+1;
                        else
                            if m>k
                                m = m-1;
                            end
                            
                            res(end+1).type   = 'artefact_flat';
                            res(end).value    = char(D.chanlabels(chanind(j)));
                            res(end).time     = D.time(onsetsamples(k)+1) - D.time(1) + D.trialonset(i);
                            res(end).duration = (onsetsamples(m) + diffs(onsets(m))-onsetsamples(k))/D.fsample;
                            
                            k = m+1;
                            m = k;
                        end
                    else
                        ind1 = onsetsamples(k);
                        ind2 = length(dat);
                        if (sum(dat(ind1:ind2))/(ind2-ind1+1))<0.5
                            res(end+1).type   = 'artefact_flat';
                            res(end).value    = char(D.chanlabels(chanind(j)));
                            res(end).time     = D.time(onsetsamples(k)+1) - D.time(1) + D.trialonset(i);
                            res(end).duration = (length(dat)-onsetsamples(k)+2)/D.fsample;
                        else                            
                            m = m-1;                            
                            
                            res(end+1).type   = 'artefact_flat';
                            res(end).value    = char(D.chanlabels(chanind(j)));
                            res(end).time     = D.time(onsetsamples(k)+1) - D.time(1) + D.trialonset(i);
                            res(end).duration = (onsetsamples(m) + diffs(onsets(m))-onsetsamples(k))/D.fsample;                                                    
                        end
                        break;
                    end
                end
            end
            if isequal(D.type, 'continuous')
                if any(Ibar == j), spm_progress_bar('Set', j); end
            end
        end
        if ~isempty(res)
            ev = D.events(i);
            if iscell(ev)
                ev = ev{1};
            end
            
            if ~S.append
                ev(strmatch('artefact_flat', {ev.type})) = [];
            end
                            
            D = events(D, i, spm_cat_struct(ev, res));
        end
        
        if ~isequal(D.type, 'continuous')
            if any(Ibar == i), spm_progress_bar('Set', i); end
        end
    end
    
    spm_progress_bar('Clear');
    
    res = D;
end

spm('FigName','M/EEG flat data detection: done');
