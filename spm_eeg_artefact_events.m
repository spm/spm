function res = spm_eeg_artefact_events(S)
% Plugin for spm_eeg_artefact for rejection based on events
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
% $Id: spm_eeg_artefact_events.m 7132 2017-07-10 16:22:58Z guillaume $


%-This part if for creating a config branch that plugs into spm_cfg_eeg_artefact
% Any parameters can be specified and they are then passed to the plugin
% when it is called.
%--------------------------------------------------------------------------
if nargin == 0
    artefacts      = cfg_const;
    artefacts.tag  = 'artefacts';
    artefacts.name = 'All artefact events';
    artefacts.val  = {1};
    artefacts.help = {''};
    
    eventlist        = cfg_files;
    eventlist.tag    = 'eventlist';
    eventlist.name   = 'Load event list';
    eventlist.filter = 'mat';
    eventlist.num    = [1 1];
    eventlist.help   = {'Select events list file.'};
    
    whatevents        = cfg_choice;
    whatevents.tag    = 'whatevents';
    whatevents.name   = 'What events to use?';
    whatevents.values = {artefacts, eventlist};
    whatevents.val    = {artefacts};
    whatevents.help   = {''};
    
    events      = cfg_branch;
    events.tag  = 'events';
    events.name = 'Reject based on events';
    events.val  = {whatevents};
    events.help = {''};
    
    res = events;
    
    return
end

SVNrev = '$Rev: 7132 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG event-based rejection');

if isequal(S.mode, 'mark')
    error('Only reject mode is supported by this plug-in.');
end

D = spm_eeg_load(S.D);

chanind  = S.chanind;
res = zeros(D.nchannels, D.ntrials);

if isfield(S.whatevents, 'eventlist')
    ev = getfield(load(char(S.whatevents.eventlist)), 'events');
else
    ev = [];
end

%-Artefact detection
%--------------------------------------------------------------------------

spm_progress_bar('Init', D.ntrials, 'Trials checked');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
else Ibar = [1:D.ntrials]; end

for i = 1:D.ntrials
    if isempty(ev)
        res(chanind, i) = squeeze(any(badsamples(D, chanind, ':', i), 2));
    else
        cev = D.events(i);
        if iscell(cev)
            cev = cev{1};
        end
        if ~isempty(cev)
            [sel1, sel2] = spm_match_str({ev.type}, {cev.type});
            for k = 1:length(sel1)
                sel3 = strmatch(ev(sel1(k)).type, {cev(sel2).type});
                if strncmp('artefact_', ev(sel1(k)).type, 9)
                    sel4 = spm_match_str(D.chanlabels(chanind), {cev(sel2(sel3)).value});
                    res(chanind(sel4), i) = 1;
                elseif isempty(ev(sel1(k)).value)
                    if any(cellfun('isempty', {cev(sel2(sel3)).value}))
                        res(chanind, i) = 1;
                    end
                elseif isnumeric(ev(sel1(k)).value)
                    sel4 = find(cellfun(@isnumeric, {cev(sel2(sel3)).value}));
                    if ~isempty(sel3) && (any(ev(sel1(k)).value == [cev(sel2(sel3(sel4))).value]))
                        res(chanind, i) = 1;
                    end
                elseif ischar(ev(sel1(k)).value)
                    sel4 = find(cellfun(@ischar, {cev(sel2(sel3)).value}));
                    if ~isempty(strmatch(ev(sel1(k)).value, {cev(sel2(sel3(sel4))).value}))
                        res(chanind, i) = 1;
                    end
                end
            end
        end
    end
    if ismember(i, Ibar), spm_progress_bar('Set', i); end
end

spm_progress_bar('Clear');

spm('FigName', 'M/EEG threshold channels: done');
