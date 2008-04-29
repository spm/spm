function [trl, conditionlabels] = spm_eeg_definetrial(S)
% Function for definition of trials based on events
% FORMAT S = spm_eeg_definetrial(S)
% S - existing configuration struct (optional)
% Fields of S:
%   S.event - event struct  (optional)
%   S.fsample - sampling rate
%   S.dataset - raw dataset (events and fsample can be read from there if absent)
%   S.timeonset - time of the first sample in the data (default - 0)
%   S.pretrig - pre-trigger time in ms
%   S.posttrig - post-trigger time in ms.
%   S.trialdef - structure array for trial definition with fields
%       S.trialdef.conditionlabel - string label for the condition
%       S.trialdef.eventtype  - string
%       S.trialdef.eventvalue  - string, numeric or empty
%   S.review - 1 - review individual trials after selection (0 - not)
% OUTPUT:
%   trl - Nx3 matrix [start end offset]
%   conditionlabels - Nx1 cell array of strings, label for each trial
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Robert Oostenveld
% $Id: spm_eeg_definetrial.m 1516 2008-04-29 15:03:26Z vladimir $

if nargin == 0
    S = [];
end

% ------------- Check inputs

Fig = spm_figure('GetWin','Interactive');
clf(Fig);

if ~isfield(S, 'event') || ~isfield(S, 'fsample')
    if ~isfield(S, 'dataset')
        S.dataset = spm_select(1, '\.*', 'Select M/EEG data file');
    end
    event = fileio_read_event(S.dataset, 'detectflank', 'both');
    hdr = fileio_read_header(S.dataset, 'fallback', 'biosig');
    S.fsample = hdr.Fs;
else
    event = S.event;
end

if ~isfield(event, 'time')
    for i = 1:numel(event)
        event(i).time = event(i).sample./S.fsample;
    end
end

if ~isfield(event, 'sample')
    if ~isfield(S, 'timeonset')
        S.timeonset = 0;
    end

    for i = 1:numel(event)
        if S.timeonset == 0
            event(i).sample = event(i).time*S.fsample;
        else
            event(i).sample = (event(i).time-S.timeonset)*S.fsample+1;
        end
    end
end

if ~isfield(S, 'pretrig')
    S.pretrig = spm_input('Start of trial in PST [ms]', '+1', 'r', '', 1);
end

if ~isfield(S, 'posttrig')
    S.posttrig = spm_input('End of trial in PST [ms]', '+1', 'r', '', 1);
end

if ~(isfield(S, 'trialdef'))
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

% ------------- Build trl based on selected events

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
        trlbeg = event(sel(j)).sample + trloff;
        trldur = round(0.001*(-S.pretrig+S.posttrig)*S.fsample);
        trlend = trlbeg + trldur;
        % add the beginsample, endsample and offset of this trial to the list
        trl = [trl; trlbeg trlend trloff];
        conditionlabels{end+1} = S.trialdef(i).conditionlabel;
    end
end

% sorts the trl in right temporal order
[junk sortind]=sort(trl(:,1));
trl=trl(sortind, :);
conditionlabels = conditionlabels(sortind);

% ------------- Review selected trials

if ~isfield(S, 'review')
    S.review = spm_input('Review individual trials?','+1','yes|no',[1 0], 0);
end

if S.review

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

% ------------- Create trial definition file

if ~isfield(S, 'save')
    S.save = spm_input('Save trial definition?','+1','yes|no',[1 0], 0);
end

if S.save
    [trlfilename, trlpathname] = uiputfile( ...
        {'*.mat', 'MATLAB File (*.mat)'}, 'Save trial definition as');

    save(fullfile(trlpathname, trlfilename), 'trl', 'conditionlabels');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that allows the user to select an event using gui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selected=select_event_ui(event)

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
