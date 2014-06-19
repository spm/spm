function res = spm_eeg_artefact_nans(S)
% Plugin for spm_eeg_artefact doing NaN detection.
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
% Copyright (C) 2011-2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_artefact_nans.m 6060 2014-06-19 13:31:19Z vladimir $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_artefact
% Any parameters can be specified and they are then passed to the plugin
% when it's called.
%--------------------------------------------------------------------------
if nargin == 0
    nans = cfg_branch;
    nans.tag = 'nans';
    nans.name = 'Detect NaNs';
    nans.val = {};
    
    res = nans;
    
    return
end

SVNrev = '$Rev: 6060 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG NaN detection');

D = spm_eeg_load(S.D);

chanind  =  S.chanind;
res = zeros(D.nchannels, D.ntrials);

if isequal(S.mode, 'reject')
    res = zeros(D.nchannels, D.ntrials);
    
    %-Artefact detection
    %--------------------------------------------------------------------------
    
    spm_progress_bar('Init', D.ntrials, 'Trials checked');
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
    else Ibar = [1:D.ntrials]; end
    
    for i = 1:D.ntrials
        for j = 1:length(chanind)
            if any(isnan(squeeze(D(chanind(j), :, i))))
                res(chanind(j), i) = 1;
            end
        end
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
            dat  = ~isnan(squeeze(D(chanind(j), :, i)));
            if  sum(dat)/length(dat)<(1-S.badchanthresh)
                res(end+1).type   = 'artefact_nan';
                res(end).value    = char(D.chanlabels(chanind(j)));
                res(end).time     = D.trialonset(i);
                res(end).duration = D.time(end) - D.time(1) + 1;
            else
                tmp  = find(dat);
                diffs = diff([0 tmp D.nsamples]);
                onsets = find(diffs>1);
                
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
                            
                            res(end+1).type   = 'artefact_nan';
                            res(end).value    = char(D.chanlabels(chanind(j)));
                            res(end).time     = D.time(onsetsamples(k)+1) - D.time(1) + D.trialonset(i);
                            res(end).duration = (onsetsamples(m) + diffs(onsets(m))-onsetsamples(k)-1)/D.fsample;
                            
                            k = m+1;
                            m = k;
                        end
                    else
                        ind1 = onsetsamples(k);
                        ind2 = length(dat);
                        if (sum(dat(ind1:ind2))/(ind2-ind1+1))<0.5
                            res(end+1).type   = 'artefact_nan';
                            res(end).value    = char(D.chanlabels(chanind(j)));
                            res(end).time     = D.time(onsetsamples(k)+1) - D.time(1) + D.trialonset(i);
                            res(end).duration = (length(dat)-onsetsamples(k)+2)/D.fsample;
                        else
                            m = m-1;
                            
                            res(end+1).type   = 'artefact_nan';
                            res(end).value    = char(D.chanlabels(chanind(j)));
                            res(end).time     = D.time(onsetsamples(k)+1) - D.time(1) + D.trialonset(i);
                            res(end).duration = (onsetsamples(m) + diffs(onsets(m))-onsetsamples(k)-1)/D.fsample;
                        end
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
                ev(strmatch('artefact_nan', {ev.type})) = [];
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

spm('FigName','M/EEG NaN detection: done');