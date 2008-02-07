function [event, recoded_events ]=spm_eeg_recode_events(event)
% This function modifies the event structure of Fieldtrip to enable SPM to
% handle it. In particular events whose value is empty are assigned the
% value Inf (for numeric) or 'Inf' (for strings). In addition the function
% assigns each compbination of type and value a unique code that is later
% used in SPM file.
%
% We would like to assign codes to events that one one hand would be
% consistent between files, but on the other hand would not clash and
% handle both string and numeric type/value cases. This is implemented
% by trying to convert the original type/value to a number with a hash
% function and then checking for the unlikely event that there is a
% previous code with the same value and assigning a new code in the
% unlikely event that happens.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging


recoded_events=[];

if isempty(event)
    warning('No events were found');
    return;
end

eventtype = unique({event.type});
Neventtype = length(eventtype);



newcode=1;
allcodes=[];

for i=1:Neventtype
    type_sel = find(strcmp(eventtype{i}, {event.type}));

    emptyval=find(cellfun('isempty', {event(type_sel).value}));

    if ~any(cellfun('isclass', {event(type_sel).value}, 'char'))

        % Incompatible with old Matlab versions
        %[event(type_sel(emptyval)).value]=deal(Inf);

        for n=1:length(emptyval)
            event(type_sel(emptyval(n))).value=Inf;
        end

        eventvalue = unique([event(type_sel).value]);
    else
        if ~isempty(strmatch('Inf', {event(type_sel).value},'exact'))
            % It's a very unlikely scenario but ...
            warning('Event value ''Inf'' cannot be handled by GUI selection. Mistakes are possible.')
        end

        % Incompatible with old Matlab versions
        % [event(type_sel(emptyval)).value]=deal('Inf');

        for n=1:length(emptyval)
            event(type_sel(emptyval(n))).value='Inf';
        end

        eventvalue = unique({event(type_sel).value});

        if ~iscell(eventvalue)
            eventvalue={eventvalue};
        end
    end
    for j=1:length(eventvalue)

        eventcode=hash(eventtype(i), eventvalue(j));

        if any(allcodes==eventcode)
            warning('Codes clushed. Assigning a new code.');
            eventcode=newcode;
            newcode=newcode+1;
        end

        allcodes=[allcodes eventcode];


        if isnumeric(eventvalue)
            value_sel=find(eventvalue(j)== [event(type_sel).value]);
        else
            value_sel=find(strcmp(eventvalue{j}, {event(type_sel).value}));
        end

        recoded_events=[recoded_events; eventcode*ones(length(value_sel) ,1) [event(type_sel(value_sel)).sample]'];

        %Incompatible with old Matlab versions
        %[event(type_sel(value_sel)).spmcode]=deal(eventcode);
        for n=1:length(value_sel)
            event(type_sel(value_sel(n))).spmcode=eventcode;
        end
    end % Loop over values of a particular type
end % Loop over types

recoded_events=sortrows(recoded_events,2);


function code=hash(varargin)
% This function converts the inputs (string or numbers) to a hopefully
% unique large number that is used as a code.

% This value is added to the codes so that they will never clash with the
% new codes. (Unless there are many new code which is very unlikely).
MAX_NEW_CODE=100;

code=0;
for n=1:nargin
    if iscell(varargin{n})
        varargin{n}=varargin{n}{1};
    end
    if ~isempty(varargin{n}) & ~any(isinf(varargin{n}))
        addcode=abs(sum(double(varargin{n})));
        % Handle some pathological cases
        if (addcode==0)
            sum(double('spm_zero'));
        end
        % In some cases the value is time. Then there are only fractional
        % differences between different events. This is done so that
        % the same code will not be assigned
        addcode=fix(addcode)+1/(0.1+addcode);
        code=10^(n-1)*code+addcode;
    end
end

code=round(code+MAX_NEW_CODE);


