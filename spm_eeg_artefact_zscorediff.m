function res = spm_eeg_artefact_zscorediff(S)
% Plugin for spm_eeg_artefact doing z-score thresholding on the diff time series
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
% Copyright (C) 2013-2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_artefact_zscorediff.m 7132 2017-07-10 16:22:58Z guillaume $


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
    threshold.val     = {5};
    threshold.help    = {'Threshold value (in stdev).'};
    
    excwin         = cfg_entry;
    excwin.tag     = 'excwin';
    excwin.name    = 'Excision window';
    excwin.strtype = 'r';
    excwin.num     = [1 1];
    excwin.val     = {100};
    excwin.help    = {'Window (in ms) to mark as bad around each event (for mark mode only), 0 - do not mark data as bad.'};
    
    zscorediff      = cfg_branch;
    zscorediff.tag  = 'zscorediff';
    zscorediff.name = 'Threshold z-scored difference data';
    zscorediff.val  = {threshold, excwin};
    zscorediff.help = {''};
    
    res = zscorediff;
    
    return
end

SVNrev = '$Rev: 7132 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG zscore diff thresholding');

D = spm_eeg_load(S.D);

chanind = S.chanind;
threshold = S.threshold;

if isequal(S.mode, 'reject') && isequal(D.type, 'continuous')
    error('Rejection mode not for continuous data.');
end

res = zeros(D.nchannels, D.ntrials);

%-Artefact detection
%--------------------------------------------------------------------------

spm_progress_bar('Init', length(chanind), 'Channels checked');
if length(chanind) > 100, Ibar = floor(linspace(1, length(chanind),100));
else Ibar = [1:length(chanind)]; end

for j = 1:length(chanind)
    dat = squeeze(D(chanind(j), :, :));
    if size(dat, 1) == 1
        dat = dat';
    end
    
    dat = diff(dat);
    
    zsdat = dat./repmat(std(dat), size(dat, 1), 1);
    
    zsdat = [zeros(1, size(zsdat, 2)); zsdat]; 
    
    if isequal(S.mode, 'reject')
        res(chanind(j), :) = any(abs(zsdat)>threshold);
    elseif isequal(S.mode, 'mark')
        bad  = abs(zsdat)>threshold;
        if ~any(bad(:))
            if isequal(D.type, 'continuous')
                if any(Ibar == j), spm_progress_bar('Set', j); end
            end
            continue;
        end
        
        if S.excwin>0
            excwin = ones(round(1e-3*S.excwin*D.fsample), 1);
            bad = ~~conv2(excwin, 1, double(bad), 'same');
        end
        
        for i = 1:D.ntrials
            res = [];
            
            onsets  = find(bad(:, i));
            offsets = find(~bad(:, i));
            onsets(find(diff(onsets)<2)+1) = [];
            
            if bad(end, i)
                offsets(end+1) = size(bad, 1)+1;
            end
            
            for k = 1:length(onsets)
                res(end+1).type   = 'artefact_zscorediff';
                res(end).value    = char(D.chanlabels(chanind(j)));
                res(end).time     = D.time(onsets(k)+1) - D.time(1) + D.trialonset(i);
                res(end).duration = (min(offsets(offsets>onsets(k)))-onsets(k))./D.fsample;
            end
            
            if ~isempty(res)
                ev = D.events(i);
                if iscell(ev)
                    ev = ev{1};
                end
                
                if ~S.append
                    sel1 = strmatch('artefact_zscorediff', {ev.type});
                    sel2 = strmatch(char(D.chanlabels(chanind(j))), {ev(sel1).value});
                    ev(sel1(sel2)) = [];
                end
                
                D = events(D, i, spm_cat_struct(ev, res));
            end
        end
    end
    
    if any(Ibar == j), spm_progress_bar('Set', j); end
end

spm_progress_bar('Clear');

if isequal(S.mode, 'mark')
    res = D;
end

spm('FigName','M/EEG zscore diff thresholding: done');
