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
% $Id: spm_eeg_artefact_jump.m 5613 2013-08-15 11:56:07Z vladimir $


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

SVNrev = '$Rev: 5613 $';

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
    
    
    spm_progress_bar('Init', D.ntrials, 'Trials checked');
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
    else Ibar = [1:D.ntrials]; end
    
    
    for i = 1:D.ntrials
        
        bad  = abs(diff(squeeze(D(chanind, :, i)), [], 2))>threshold;
        
        if ~any(bad(:))
            if ismember(i, Ibar), spm_progress_bar('Set', i); end
            continue;
        end
        
        if S.excwin>0
            excwin = ones(round(5e-4*S.excwin*D.fsample), 1);
            bad = ~~conv2(excwin, 1, double(bad), 'same');
        end
        
        res = [];
        for j = 1:length(chanind)
            onsets  = find(bad(j, :));
            offsets = find(~bad(j, :));
            onsets(find(diff(onsets)<2)+1) = [];
            
            for k = 1:length(onsets)
                res(end+1).type   = 'artefact_jump';
                res(end).value    = char(D.chanlabels(chanind(j)));
                res(end).time     = D.time(onsets(k)) - D.time(1) + D.trialonset(i);
                res(end).duration = (min(offsets(offsets>onsets(k)))-onsets(k)+1)./D.fsample;
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
        
        
        if ismember(i, Ibar), spm_progress_bar('Set', i); end
    end
    
    spm_progress_bar('Clear');
    
    res = D;
end

spm('FigName','M/EEG jump detection: done');