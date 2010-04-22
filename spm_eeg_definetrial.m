function [trl, conditionlabels, S] = spm_eeg_definetrial(S)
% Definition of trials based on events
% FORMAT[trl, conditionlabels, S]  = spm_eeg_definetrial(S)
% S                - input structure (optional)
% (optional) fields of S:
%   S.event        - event struct  (optional)
%   S.fsample      - sampling rate
%   S.dataset      - raw dataset (events and fsample can be read from there if absent)
%   S.inputformat  - data type (optional) to force the use of specific data reader
%   S.timeonset    - time of the first sample in the data [default: 0]
%   S.pretrig      - pre-trigger time in ms
%   S.posttrig     - post-trigger time in ms
%   S.trialdef     - structure array for trial definition with fields (optional)
%       S.trialdef.conditionlabel - string label for the condition
%       S.trialdef.eventtype      - string
%       S.trialdef.eventvalue     - string, numeric or empty
%   S.reviewtrials - review individual trials after selection (yes/no: 1/0)
%   S.save         - save trial definition (yes/no: 1/0)
% OUTPUT:
%   trl            - Nx3 matrix [start end offset]
%   conditionlabels - Nx1 cell array of strings, label for each trial
%   S              - modified configuration structure (for history)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Robert Oostenveld
% $Id: spm_eeg_definetrial.m 3833 2010-04-22 14:49:48Z vladimir $


SVNrev = '$Rev: 3833 $';

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','M/EEG trial definition');

%-Get input parameters
%--------------------------------------------------------------------------
try
    S.inputformat;
catch
    S.inputformat = [];
end

if ~isfield(S, 'event') || ~isfield(S, 'fsample')
    if ~isfield(S, 'dataset')
        S.dataset = spm_select(1, '\.*', 'Select M/EEG data file');
    end
    
    hdr = ft_read_header(S.dataset, 'fallback', 'biosig', 'headerformat', S.inputformat);
    S.fsample = hdr.Fs;
    
    event = ft_read_event(S.dataset, 'detectflank', 'both', 'eventformat', S.inputformat);

    if ~isempty(strmatch('UPPT001', hdr.label))
        % This is s somewhat ugly fix to the specific problem with event
        % coding in FIL CTF. It can also be useful for other CTF systems where the
        % pulses in the event channel go downwards.
        fil_ctf_events = ft_read_event(S.dataset, 'detectflank', 'down', 'type', 'UPPT001', 'trigshift', -1, 'eventformat', S.inputformat);
        if ~isempty(fil_ctf_events)
            [fil_ctf_events(:).type] = deal('FIL_UPPT001_down');
            event = cat(1, event(:), fil_ctf_events(:));
        end
    end
    
    
    if ~isempty(strmatch('UPPT002', hdr.label))
        % This is s somewhat ugly fix to the specific problem with event
        % coding in FIL CTF. It can also be useful for other CTF systems where the
        % pulses in the event channel go downwards.
        fil_ctf_events = ft_read_event(S.dataset, 'detectflank', 'down', 'type', 'UPPT002', 'trigshift', -1, 'eventformat', S.inputformat);
        if ~isempty(fil_ctf_events)
            [fil_ctf_events(:).type] = deal('FIL_UPPT002_down');
            event = cat(1, event(:), fil_ctf_events(:));
        end
    end


    % This is another FIL-specific fix that will hopefully not affect other sites
    if isfield(hdr, 'orig') && isfield(hdr.orig, 'VERSION') && isequal(uint8(hdr.orig.VERSION),uint8([255 'BIOSEMI']))
        ind = strcmp('STATUS', {event(:).type});
        val = [event(ind).value];
        if any(val>255)
            bytes  = dec2bin(val);
            bytes  = bytes(:, end-7:end);
            bytes  = flipdim(bytes, 2);
            val    = num2cell(bin2dec(bytes));
            [event(ind).value] = deal(val{:});
        end
    end

else
    event = S.event;
end

if ~isfield(S, 'timeonset')
    S.timeonset = 0;
end

if ~isfield(event, 'time')
    for i = 1:numel(event)
        if S.timeonset == 0
            event(i).time = event(i).sample./S.fsample;
        else
            event(i).time = (event(i).sample - 1)./S.fsample + S.timeonset;
        end
    end
end

if ~isfield(event, 'sample')
    for i = 1:numel(event)
        if S.timeonset == 0
            event(i).sample = event(i).time*S.fsample;
        else
            event(i).sample = (event(i).time-S.timeonset)*S.fsample+1;
        end
        
        event(i).sample = round(event(i).sample);
    end
end

if isempty(event)
    error('No event information was found in the input');
end

if ~isfield(S, 'pretrig')
    S.pretrig = spm_input('Start of trial in PST [ms]', '+1', 'r', '', 1);
end

if ~isfield(S, 'posttrig')
    S.posttrig = spm_input('End of trial in PST [ms]', '+1', 'r', '', 1);
end

if ~isfield(S, 'trialdef')
    S.trialdef = [];
    ncond = spm_input('How many conditions?', '+1', 'n', '1');
    for i = 1:ncond
        OK = 0;
        pos = '+1';
        while ~OK
            conditionlabel = spm_input(['Label of condition ' num2str(i)], pos, 's');
            selected = select_event_ui(event);
            if isempty(conditionlabel) || isempty(selected)
                pos = '-1';
            else
                for j = 1:size(selected, 1)
                    S.trialdef = [S.trialdef ...
                        struct('conditionlabel', conditionlabel, ...
                        'eventtype', selected{j, 1}, ...
                        'eventvalue', selected{j, 2})];
                    OK=1;
                end
            end
        end
    end
end

for i = 1:length(S.trialdef)
    if ~isfield(S.trialdef(i),'trlshift')
        trlshift(i) = 0;
    else
        trlshift(i) = round(S.trialdef(i).trlshift * S.fsample/1000); % assume passed as ms
    end
end

%-Build trl based on selected events
%--------------------------------------------------------------------------
trl = [];
conditionlabels = {};
for i=1:numel(S.trialdef)

    if ischar(S.trialdef(i).eventvalue)
        % convert single string into cell-array, otherwise the intersection does not work as intended
        S.trialdef(i).eventvalue = {S.trialdef(i).eventvalue};
    end

    sel = [];
    % select all events of the specified type and with the specified value
    for j=find(strcmp(S.trialdef(i).eventtype, {event.type}))
        if isempty(S.trialdef(i).eventvalue)
            sel = [sel j];
        elseif ~isempty(intersect(event(j).value, S.trialdef(i).eventvalue))
            sel = [sel j];
        end
    end

    for j=1:length(sel)
        % override the offset of the event
        trloff = round(0.001*S.pretrig*S.fsample);        
        % also shift the begin sample with the specified amount
        if ismember(event(sel(j)).type, {'trial', 'average'})
            % In case of trial events treat the 0 time point as time of the
            % event rather than the beginning of the trial 
            trlbeg = event(sel(j)).sample - event(sel(j)).offset + trloff;
        else
            trlbeg = event(sel(j)).sample + trloff;
        end
        trldur = round(0.001*(-S.pretrig+S.posttrig)*S.fsample);
        trlend = trlbeg + trldur;
        
        % Added by Rik in case wish to shift triggers (e.g, due to a delay
        % between trigger and visual/auditory stimulus reaching subject).
        trlbeg = trlbeg + trlshift(i);
        trlend = trlend + trlshift(i);
        
        % add the beginsample, endsample and offset of this trial to the list
        trl = [trl; trlbeg trlend trloff];
        conditionlabels{end+1} = S.trialdef(i).conditionlabel;
    end
end

%-Sort the trl in right temporal order
%--------------------------------------------------------------------------
[junk, sortind] = sort(trl(:,1));
trl             = trl(sortind, :);
conditionlabels = conditionlabels(sortind);

%-Review selected trials
%--------------------------------------------------------------------------
if ~isfield(S, 'reviewtrials')
    S.reviewtrials = spm_input('Review individual trials?','+1','yes|no',[1 0], 0);
end

if S.reviewtrials

    eventstrings=cell(size(trl,1),1);
    for i=1:size(trl,1)
        eventstrings{i}=[num2str(i) ' Label: ' conditionlabels{i} ' Time (sec): ' num2str((trl(i, 1)- trl(i, 3))./S.fsample)];
    end

    selected = find(trl(:,1)>0);

    [indx OK] = listdlg('ListString', eventstrings, 'SelectionMode', 'multiple', 'InitialValue', ...
        selected, 'Name', 'Select events', 'ListSize', [300 300]);

    if OK
        trl=trl(indx, :);
    end
end

%-Create trial definition file
%--------------------------------------------------------------------------
if ~isfield(S, 'save')
    S.save = spm_input('Save trial definition?','+1','yes|no',[1 0], 0);
end

if S.save
    [trlfilename, trlpathname] = uiputfile( ...
        {'*.mat', 'MATLAB File (*.mat)'}, 'Save trial definition as');

    save(fullfile(trlpathname, trlfilename), 'trl', 'conditionlabels');
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG trial definition: done');

%==========================================================================
% select_event_ui
%==========================================================================
function selected = select_event_ui(event)
% Allow the user to select an event using GUI

selected={};

if isempty(event)
    fprintf('no events were found\n');
    return
end

eventtype = unique({event.type});
Neventtype = length(eventtype);

% Two lists are built in parallel
settings={}; % The list of actual values to be used later
strsettings={}; % The list of strings to show in the GUI

for i=1:Neventtype
    sel = find(strcmp(eventtype{i}, {event.type}));

    numind = find(...
        cellfun('isclass', {event(sel).value}, 'double') & ...
        ~cellfun('isempty', {event(sel).value}));

    charind = find(cellfun('isclass', {event(sel).value}, 'char'));

    emptyind = find(cellfun('isempty', {event(sel).value}));

    if ~isempty(numind)
        numvalue = unique([event(sel(numind)).value]);
        for j=1:length(numvalue)
            ninstances = sum([event(sel(numind)).value] == numvalue(j));
            strsettings=[strsettings; {['Type: ' eventtype{i} ' ; Value: ' num2str(numvalue(j)) ...
                ' ; ' num2str(ninstances) ' instances']}];
            settings=[settings; [eventtype(i), {numvalue(j)}]];
        end
    end

    if ~isempty(charind)
        charvalue = unique({event(sel(charind)).value});
        if ~iscell(charvalue)
            charvalue = {charvalue};
        end
        for j=1:length(charvalue)
            ninstances = length(strmatch(charvalue{j}, {event(sel(charind)).value}));
            strsettings=[strsettings; {['Type: ' eventtype{i} ' ; Value: ' charvalue{j}...
                ' ; ' num2str(ninstances) ' instances']}];
            settings=[settings; [eventtype(i), charvalue(j)]];
        end
    end

    if ~isempty(emptyind)
        strsettings=[strsettings; {['Type: ' eventtype{i} ' ; Value: ; ' ...
            num2str(length(emptyind)) ' instances']}];
        settings=[settings; [eventtype(i), {[]}]];
    end
end

[selection ok]= listdlg('ListString',strsettings, 'SelectionMode', 'multiple', 'Name', 'Select event', 'ListSize', [400 300]);

if ok
    selected=settings(selection, :);
else
    selected={};
end
