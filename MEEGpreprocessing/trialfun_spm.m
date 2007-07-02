function [trl, event] = trialfun_spm(cfg);

% TRIALFUN_SPM  determines trials/segments in the data that are
% interesting for analysis, using the general event structure returned
% by read_fcdc_event. This function is independent of the dataformat
%
% The trialdef structure can contain the following specifications
%   cfg.trialdef.eventtype  = 'string'
%   cfg.trialdef.eventvalue = number, string or list with numbers or strings
%   cfg.trialdef.prestim    = latency in seconds (optional)
%   cfg.trialdef.poststim   = latency in seconds (optional)
%
% If you want to read all data from a continous file in segments, you can specify
%    cfg.trialdef.triallength = duration in seconds (can be Inf)
%    cfg.trialdef.ntrials     = number of trials
%
% If you specify
%   cfg.trialdef.eventtype  = '?'
% a list with the events in your datafile will be displayed on screen.

% Copyright (C) 2005-2007, Robert Oostenveld
%
% $Log: trialfun_general.m,v $
% Revision 1.6  2007/02/07 16:30:12  roboos
% do something with trials that overlap with teh beginning of the file (in case of gui)
%
% Revision 1.5  2007/02/07 07:30:30  roboos
% added gui option, thanks to Vladimir
%
% Revision 1.4  2007/01/09 09:47:25  roboos
% added warning about this function being deprecated, includes link to website documentation
%

% some events do not require the specification a type, pre or poststim period
% in that case it is more convenient not to have them, instead of making them empty


if isfield(cfg.trialdef, 'eventvalue')  && isempty(cfg.trialdef.eventvalue   ), cfg.trialdef = rmfield(cfg.trialdef, 'eventvalue' ); end
if isfield(cfg.trialdef, 'prestim')     && isempty(cfg.trialdef.prestim      ), cfg.trialdef = rmfield(cfg.trialdef, 'prestim'    ); end
if isfield(cfg.trialdef, 'poststim')    && isempty(cfg.trialdef.poststim     ), cfg.trialdef = rmfield(cfg.trialdef, 'poststim'   ); end
if isfield(cfg.trialdef, 'triallength') && isempty(cfg.trialdef.triallength  ), cfg.trialdef = rmfield(cfg.trialdef, 'triallength'); end
if isfield(cfg.trialdef, 'ntrials')     && isempty(cfg.trialdef.ntrials      ), cfg.trialdef = rmfield(cfg.trialdef, 'ntrials'    ); end

if isfield(cfg.trialdef, 'triallength')
    % reading all segments from a continuous fie is incompatible with any other option
    try, cfg.trialdef = rmfield(cfg.trialdef, 'eventvalue'); end
    try, cfg.trialdef = rmfield(cfg.trialdef, 'prestim'   ); end
    try, cfg.trialdef = rmfield(cfg.trialdef, 'poststim'  ); end
    if ~isfield(cfg.trialdef, 'ntrials')
        if isinf(cfg.trialdef.triallength)
            cfg.trialdef.ntrials = 1;
        else
            cfg.trialdef.ntrials = inf;
        end
    end
end

% read the header and event information
hdr = read_header(cfg.dataset);
try
    event = read_event(cfg.dataset);
catch
    event = [];
end

event=spm_eeg_recode_events(event);

% for the following, the trials do not depend on the events in the data
if isfield(cfg.trialdef, 'triallength')
    if isinf(cfg.trialdef.triallength)
        % make one long trial with the complete continuous data in it
        trl = [1 hdr.nSamples*hdr.nTrials 0];
    elseif isinf(cfg.trialdef.ntrials)
        % cut the continous data into as many segments as possible
        nsamples = round(cfg.trialdef.triallength*hdr.Fs);
        trlbeg   = 1:nsamples:(hdr.nSamples*hdr.nTrials - nsamples + 1);
        trlend   = trlbeg + nsamples - 1;
        offset   = zeros(size(trlbeg));
        trl = [trlbeg(:) trlend(:) offset(:)];
    else
        % make the pre-specified number of trials
        nsamples = round(cfg.trialdef.triallength*hdr.Fs);
        trlbeg   = (0:(cfg.trialdef.ntrials-1))*nsamples + 1;
        trlend   = trlbeg + nsamples - 1;
        offset   = zeros(size(trlbeg));
        trl = [trlbeg(:) trlend(:) offset(:)];
    end
    return
end

if strcmp(cfg.trialdef.eventtype, '?')
    % no trials should be added, show event information using subfunction and exit
    show_event(event);
    return
end

if strcmp(cfg.trialdef.eventtype, 'gui') || (isfield(cfg.trialdef, 'eventvalue') && strcmp(cfg.trialdef.eventvalue, 'gui'))
    selected = select_event(event, cfg.trialdef);
    usegui = 1;
else
    selected={cfg.trialdef.eventtype, cfg.trialdef.eventvalue};
    usegui = 0;
end

trl = [];
for s=1:size(selected, 1)

    cfg.trialdef.eventtype=selected{s,1};
    cfg.trialdef.eventvalue=selected{s,2};

    if ~isfield(cfg.trialdef, 'eventvalue')
        cfg.trialdef.eventvalue = [];
    elseif ischar(cfg.trialdef.eventvalue)
        % convert single string into cell-array, otherwise the intersection does not work as intended
        cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
    end

    sel = [];
    % select all events of the specified type and with the specified value
    for i=find(strcmp(cfg.trialdef.eventtype, {event.type}))
        if isempty(cfg.trialdef.eventvalue)
            sel = [sel i];
        elseif ~isempty(intersect(event(i).value, cfg.trialdef.eventvalue))
            sel = [sel i];
        end
    end

    if usegui
        % Checks whether offset and duration are defined for all the selected
        % events and/or prestim/poststim are defined in trialdef.
        if (any(cellfun('isempty', {event(sel).offset})) || ...
                any(cellfun('isempty', {event(sel).duration}))) && ...
                ~(isfield(cfg.trialdef, 'prestim') && isfield(cfg.trialdef, 'poststim'))

            % If at least some of offset/duration values and prestim/poststim
            % values are missing tries to ask the user for prestim/poststim
            answer = inputdlg({'Prestimulus latency (sec)','Poststimulus latency (sec)'}, 'Enter borders');
            if isempty(answer) || any(cellfun('isempty', answer))
                error('The information in the data and cfg is insufficient to define trials.');
            else
                cfg.trialdef.prestim=str2double(answer{1});
                cfg.trialdef.poststim=str2double(answer{2});
                if isnan(cfg.trialdef.prestim) || isnan(cfg.trialdef.poststim)
                    error('Illegal input for trial borders');
                end
            end
        end % if specification is not complete
    end % if usegui

    for i=sel
        % catch empty fields in the event table and interpret them meaningfully
        if isfield(event, 'offset') && ~isempty(event(i).offset)
            % time axis has no offset relative to the event
            event(i).offset = 0;
        end

        if ~isfield(event, 'duration') || isempty(event(i).duration)
            event(i).duration = 0;
        end
        % If there is duration information in the event and the user did
        % not override the default behavior this information is used
        % instead of user defined trial borders
        if event(i).duration>1 && strcmpi(cfg.trialdef.ignorelim, 'no')
            trloff = event(i).offset;
            trlbeg = event(i).sample;
            trldur = max(event(i).duration - 1, 0);
        else
            % override the offset of the event
            trloff = round(-cfg.trialdef.prestim*hdr.Fs);
            % also shift the begin sample with the specified amount
            trlbeg = event(i).sample + trloff;
            trldur = round((cfg.trialdef.poststim+cfg.trialdef.prestim)*hdr.Fs) - 1;
        end
        trlend = trlbeg + trldur;
        % add the beginsample, endsample and offset of this trial to the list
        trl = [trl; [trlbeg trlend trloff event(i).spmcode]];
    end

end % for loop over different selected events;

% sorts the trl in right temporal order
[junk sortind]=sort(trl(:,1));
trl=trl(sortind, :);

if usegui

    % This complicated line just computes the trigger times in seconds and
    % converts them to a cell array of strings to use in the GUI
    % It is not compatible with old Matlab versions so will be replaced
    % with simpler code for the time being
    %eventstrings = cellfun(@num2str, mat2cell((trl(:, 1)- trl(:, 3))./hdr.Fs , ones(1, size(trl, 1))), 'UniformOutput', 0);

    triggertimes=(trl(:, 1)- trl(:, 3))./hdr.Fs;
    eventstrings=cell(length(triggertimes),1);
    for t=1:length(triggertimes)
        eventstrings{t}=['SPM code: ' num2str(trl(t, 4)) ' Time (sec): ' num2str(triggertimes(t))];
    end


    % Let us start with handling at least the completely unsegmented case
    % semi-automatically. The more complicated cases are better left
    % to the user.
    if hdr.nTrials==1
        selected = find(trl(:,1)>0 & trl(:,2)<=hdr.nSamples);
    else
        selected = find(trl(:,1)>0);
    end

    [indx ok] = listdlg('ListString', eventstrings, 'SelectionMode', 'multiple', 'InitialValue', ...
        selected, 'Name', 'Select events', 'ListSize', [300 300]);

    if ok
        trl=trl(indx, :);
    end
end

% Saving the trl matrix into a file
if strcmpi(cfg.trialdef.savetrl, 'yes')
    [trlfilename, trlpathname, trlfilterindex] = uiputfile( ...
        {'*.mat', 'MATLAB File (*.mat)'; '*.*',  'All Files (ASCII) (*.*)'}, ...
        'Save TRL matrix as');
    switch trlfilterindex
        case 1
            save(fullfile(trlpathname, trlfilename), 'trl');
        case 2
            save(fullfile(trlpathname, trlfilename), 'trl', '-ascii');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that shows event table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_event(event);
if isempty(event)
    fprintf('no events were found in the datafile\n');
    return
end
eventtype = unique({event.type});
Neventtype = length(eventtype);
if Neventtype==0
    fprintf('no events were found in the datafile\n');
else
    fprintf('the following events were found in the datafile\n');
    for i=1:Neventtype
        sel = find(strcmp(eventtype{i}, {event.type}));
        try
            eventvalue = unique({event(sel).value});            % cell-array with string value
            eventvalue = sprintf('''%s'' ', eventvalue{:});     % translate into a single string
        catch
            eventvalue = unique(cell2mat({event(sel).value}));  % array with numeric values or empty
            eventvalue = num2str(eventvalue);                   % translate into a single string
        end
        fprintf('event type: ''%s'' ', eventtype{i});
        fprintf('with event values: %s', eventvalue);
        fprintf('\n');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that allows the user to select an event using gui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selected=select_event(event, trialdef)

selected={};

if isempty(event)
    fprintf('no events were found in the datafile\n');
    return
end

if strcmp(trialdef.eventtype, 'gui')
    eventtype = unique({event.type});
else
    eventtype ={trialdef.eventtype};
end

Neventtype = length(eventtype);

if Neventtype==0
    fprintf('no events were found in the datafile\n');
else
    % Two lists are built in parallel
    settings={}; % The list of actual values to be used later
    strsettings={}; % The list of strings to show in the GUI

    for i=1:Neventtype
        sel = find(strcmp(eventtype{i}, {event.type}));

        if ~any(cellfun('isclass', {event(sel).value}, 'char'))
            eventvalue = unique([event(sel).value]);
        else
            eventvalue = unique({event(sel).value});
            if ~iscell(eventvalue)
                eventvalue={eventvalue};
            end
        end
        eventcode = unique([event(sel).spmcode]);
        for j=1:length(eventvalue)

            if strcmp(trialdef.eventtype, 'gui') || strcmp(eventtype{i}, trialdef.eventtype)
                if iscell(eventvalue) && ischar(eventvalue{j}) && ~strcmp(eventvalue{j}, 'Inf')
                    settings=[settings; [eventtype(i), eventvalue(j)]];
                elseif isnumeric(eventvalue(j)) && eventvalue(j)~=Inf
                    settings=[settings; [eventtype(i), {eventvalue(j)}]];
                else
                    settings=[settings; [eventtype(i), {[]}]];
                end

                if isa(eventvalue, 'numeric')
                    strsettings=[strsettings; {['SPM code: ' num2str(eventcode(j)) ' Type: ' eventtype{i} ' ; Value: ' num2str(eventvalue(j))]}];
                else
                    strsettings=[strsettings; {['SPM code: ' num2str(eventcode(j)) ' Type: ' eventtype{i} ' ; Value: ' eventvalue{j}]}];
                end
            end
        end
    end

    if isempty(strsettings)
        fprintf('no events of the selected type were found in the datafile\n');
        return
    end

    [selection ok]= listdlg('ListString',strsettings, 'SelectionMode', 'multiple', 'Name', 'Select event', 'ListSize', [400 300]);

    if ok
        selected=settings(selection, :);
    else
        selected={};
    end
end
