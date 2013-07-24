function selected = spm_eeg_select_event_ui(event)
% Allow the user to select an event using GUI
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_select_event_ui.m 5592 2013-07-24 16:25:55Z vladimir $


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
            ninstances = length(strmatch(charvalue{j}, {event(sel(charind)).value}, 'exact'));
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

[selection,ok]= listdlg('ListString',strsettings, 'SelectionMode', 'multiple', 'Name', 'Select event', 'ListSize', [400 300]);

if ok
    selected=settings(selection, :);
else
    selected={};
end
