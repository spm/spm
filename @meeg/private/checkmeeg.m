function this = checkmeeg(this)
% Check the internal structure of meeg objects
% FORMAT this = checkmeeg(this)
% this  - the struct to check (is returned modified if necessary)
% _________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: checkmeeg.m 6817 2016-06-20 17:10:50Z vladimir $


%-Initialise data dimensions
%--------------------------------------------------------------------------
if ~isfield(this, 'Nsamples')
    this.Nsamples = 0;
end
Nsamples = this.Nsamples;

if ~isfield(this, 'channels')
    this.channels = struct([]);
end
Nchannels = length(this.channels);

if ~isfield(this, 'trials')
    this.trials = struct([]);
end
Ntrials = length(this.trials);

if ~isfield(this, 'transform')
    this.transform.ID = 'time';
end

isTF =  strncmp(this.transform.ID, 'TF', 2);

if isTF
    if ~isfield(this.transform, 'frequencies')
        error('Frequency axis must be defined for spectral dataset.');
    end
    Nfrequencies = length(this.transform.frequencies);
else
    this.transform = struct('ID', 'time');
    Nfrequencies = 0;
end

if isTF
    expected_size = [Nchannels Nfrequencies Nsamples];
else
    expected_size = [Nchannels Nsamples];
end

if Ntrials > 1
    expected_size = [expected_size Ntrials];
end

is_empty = ~all(expected_size);


%-Link to the data file if necessary
%--------------------------------------------------------------------------
if ~isfield(this, 'path') || ~exist(this.path, 'dir')
    this.path = pwd;
end

if ~isfield(this, 'data')
    this.data = [];
end

% Conversion from SPM8 format
if isfield(this.data, 'y')
    this.data = this.data.y;
end

is_linked = 0;
if isa(this.data, 'file_array')
    fname = this.data.fname;
    if ~ispc, fname = strrep(fname, '\', filesep); end
    
    [p, f, x] = fileparts(fname);
    
    % this.path takes precedence over the path in file_array to support
    % both the case when the dataset is copied (more common) and the
    % case when header saved in a mat file is used to link to a
    % datafile saved somewhere else
    fname = fullfile(this.path, [f x]);
    if ~exist(fname, 'file')
        fname = fullfile(p, [f x]);
        if ~exist(fname, 'file')
            error('Could not find the data file at either:\n  %s\nor\n  %s',...
                fullfile(this.path, [f x]),fname);
        end
    end
    
    this.data.fname = fname;
    try
        % Try reading data, i.e. check if it's a "good" filearray
        this.data(1, 1, 1);
    catch
        if ~is_empty
            scl_slope = this.data.scl_slope;
            offset    = this.data.offset;
            this.data = file_array(fname, expected_size, this.data.dtype);
            
            if  size(this.data, 1) == length(scl_slope)
                this.data.scl_slope = scl_slope;
            end
            this.data.offset = offset;
        end
    end
    
    try % once more
        this.data(1, 1, 1);
        is_linked = 1;
    catch
        warning_flexible('SPM:checkmeeg', 'Failed to link to the data file. Unlinking.');
        this.data = [];
        is_linked = 0;
    end       
end

actual_size = size(this.data); 
actual_size = [actual_size ones(1, isTF + 3 - length(actual_size))];

if ~isfield(this, 'Fsample')
    this.Fsample = 0.01;
end

if ~isfield(this, 'timeOnset')
    this.timeOnset = 0;
elseif isempty(this.timeOnset)
    if ~is_empty
        this.timeOnset = 1/this.Fsample;
    else
        this.timeOnset = 0;
    end
end

%-Check channel description
%--------------------------------------------------------------------------
if ~isfield(this, 'channels') 
    this.channels = struct([]);
end

if is_linked && ~isempty(this.channels) && (numel(this.channels) ~= actual_size(1))
    error('Channel description does not match the data.');
end
    
if is_linked
    Nchannels = actual_size(1);
else
    Nchannels = numel(this.channels);
end

if Nchannels > 0
    if ~isfield(this.channels, 'label')
        for i = 1:Nchannels
            this.channels(i).label = ['Ch' num2str(i)];
        end
    end
    
    if ~isfield(this.channels, 'bad')
        [this.channels.bad] = deal(0);
    else
        [this.channels(cellfun('isempty', {this.channels.bad})).bad] = deal(0);
    end
    
    for i = 1:Nchannels
        this.channels(i).bad = double(~~this.channels(i).bad);
    end
    
    if ~isfield(this.channels, 'type')
        [this.channels.type] = deal('Other');
    end
    
    if ~isfield(this.channels, 'X_plot2D')
        [this.channels.X_plot2D] = deal([]);
        [this.channels.Y_plot2D] = deal([]);
    end
    
    if ~isfield(this.channels, 'units')
        [this.channels.units] = deal('unknown');
    else
        [this.channels(cellfun('isempty', {this.channels.units})).units] = deal('unknown');
    end
    
    for i = 1:Nchannels
        if ~(length(this.channels(i).X_plot2D)==1 && isfinite(this.channels(i).X_plot2D))
            this.channels(i).X_plot2D = [];
        end
        
        if ~(length(this.channels(i).Y_plot2D)==1 && isfinite(this.channels(i).Y_plot2D))
            this.channels(i).Y_plot2D = [];
        end
    end
end

%-Check trial description
%--------------------------------------------------------------------------
if ~isfield(this, 'trials') 
    this.trials = struct([]);
end

if is_linked && ~isempty(this.trials) && (numel(this.trials) ~= actual_size(end))
    error('Trial description does not match the data.');
end
    
if is_linked
    Ntrials = actual_size(end);
else
    Ntrials = numel(this.trials);
end

if Ntrials > 0 
    if ~isfield(this.trials, 'label')       
        [this.trials.label] = deal('Undefined');
    elseif any(cellfun('isempty', {this.trials(:).label}))
        warning_flexible('SPM:checkmeeg', 'Some trial labels empty, assigning default.');
        [this.trials(cellfun('isempty', {this.trials(:).label})).label] = deal('Undefined');
    end
    if ~isfield(this.trials, 'bad')
        [this.trials.bad] = deal(0);
    end
    
    if ~isfield(this.trials, 'tag')
        [this.trials.tag] = deal([]);
    end
    
    if ~isfield(this.trials, 'events')
        [this.trials.events] = deal([]);
    end
    
    for i = 1:Ntrials
        
        label = this.trials(i).label;
        
        if iscell(label) && numel(label) == 1
            label = label{1};
        end
        
        if isnumeric(label)
            label = num2str(label);
        end
        
        if isa(label, 'char')
            this.trials(i).label = deblank(label);
        else
            this.trials(i).label = 'Unknown';
            warning_flexible('SPM:checkmeeg', 'Some trial labels were not strings, changing back to ''Unknown''.');
        end
        
        if  length(this.trials(i).bad)>1 || ~(this.trials(i).bad == 0 || this.trials(i).bad == 1)
            warning_flexible('SPM:checkmeeg', ['Illegal value for bad flag in trial ' num2str(i) ', resetting to zero.']);
            this.trials(i).bad = 0;
        end
        
        event = this.trials(i).events;
        
        if ~isempty(event) && ~(numel(event) == 1 && isequal(event.type, 'no events'))
            % make sure that all required elements are present
            if ~isfield(event, 'type'),     error('type field not defined for each event.');  end
            if ~isfield(event, 'time'),     error('time field not defined for each event.');  end
            if ~isfield(event, 'value'),    [event.value]    = deal([]);                     end
            if ~isfield(event, 'offset'),   [event.offset]   = deal(0);                      end
            if ~isfield(event, 'duration'), [event.duration] = deal([]);                     end
            
            
            % make sure that all numeric values are double
            for j = 1:length(event)
                if isnumeric(event(j).value)
                    event(j).value = double(event(j).value);
                end
                event(j).time      = double(event(j).time);
                event(j).offset    = double(event(j).offset);
                event(j).duration  = double(event(j).duration);
            end
            
            if ~isempty(event)
                % sort the events on the sample on which they occur
                % this has the side effect that events without time are discarded
                [junk, indx] = sort([event.time]);
                event = event(indx);
            end
        else
            event = [];
        end
        
        this.trials(i).events = event;
    end
    
    if Ntrials == 1
        this.trials.onset = this.timeOnset;
    end
    
    if ~isfield(this.trials, 'onset')
        [this.trials.onset] = deal(0);
    else
        [this.trials(cellfun('isempty', {this.trials.onset})).onset] = deal(0);
    end
    if ~isfield(this.trials, 'repl')
        [this.trials.repl] = deal(1);
    end
end

%-Check frequency axis
%--------------------------------------------------------------------------
if isTF
    if is_linked && (length(this.transform.frequencies) ~= actual_size(2))
        error('Frequency axis does not match the data.');
    end
    
    df = diff(this.transform.frequencies);
    % To avoid small numerical errors
    if length(unique(df)) > 1 && (max(diff(df))/mean(df))<0.1
        this.transform.frequencies = (0:(Nfrequencies-1))*round(100*mean(df))/100+...
            round(100*this.transform.frequencies(1))/100;
    end
end


%-Check data type
%--------------------------------------------------------------------------
if ~isfield(this, 'type') ||...
        (strcmp(this.type, 'continuous') && Ntrials>1) ||...
        (strcmp(this.type, 'evoked') && (numel(unique({this.trials.label})) ~= Ntrials))
    disp('Data type is missing or incorrect, assigning default.');
    % rule of thumb - 10 sec
    if Nsamples == 0
        this.type = 'continuous';
    elseif Ntrials==1 && (Nsamples/this.Fsample) > 10 
        this.type = 'continuous';
    elseif numel(unique({this.trials.label})) == Ntrials
        this.type = 'evoked';
    else
        this.type = 'single';
    end
end

%-Check file name
%--------------------------------------------------------------------------
if ~isfield(this, 'fname')
    if is_linked
        this.fname = [f '.mat'];
    else
        this.fname = 'spm_meeg.mat';
    end
end


%-Check sensor description
%--------------------------------------------------------------------------
if ~isfield(this, 'sensors')
    this.sensors = struct([]);
else
    if isfield(this.sensors, 'eeg')
        if isempty(this.sensors.eeg)
            this.sensors = rmfield(this.sensors, 'eeg');
        else
            try
                % This can be removed in the future
                if ~isempty(strmatch('uV', this.sensors.eeg.chanunit))
                    this.sensors.eeg = rmfield(this.sensors.eeg, 'chanunit');
                end
            end
            this.sensors.eeg = ft_datatype_sens(this.sensors.eeg, 'amplitude', 'V', 'distance', 'mm');        
        end
    end
    if isfield(this.sensors, 'meg')
        if isempty(this.sensors.meg)           
            this.sensors = rmfield(this.sensors, 'meg');
        else
            try
                % This can be removed in the future
                if ~isempty(strmatch('fT', this.sensors.meg.chanunit))
                    this.sensors.meg = rmfield(this.sensors.meg, 'chanunit');
                end
            end
            this.sensors.meg = ft_datatype_sens(this.sensors.meg, 'amplitude', 'T', 'distance', 'mm');
        end
    end
end

%-Check other fields
%--------------------------------------------------------------------------
if ~isfield(this, 'fiducials')
   this.fiducials = struct([]);
else
   this.fiducials  = ft_struct2double(fixpnt(this.fiducials));
end

if ~isfield(this, 'artifacts')
    this.artifacts = struct([]);
end

if ~isfield(this, 'other')
    this.other = struct([]);
end

if ~isfield(this, 'condlist')
    if isfield(this.other, 'condlist')
        this.condlist = this.other.condlist;
        this.other = rmfield(this.other, 'condlist');
    else
        this.condlist = {};
    end
end

if ~isempty(this.condlist)
    [junk, ind] = unique(this.condlist, 'first');
    this.condlist =  this.condlist(sort(ind));
end

if ~isfield(this, 'history')
    this.history = struct([]);
end

if ~isfield(this, 'montage') ||  ~isfield(this.montage,'M')
    this.montage = struct('M',[],'Mind',0);
else    
    if this.montage.Mind > numel(this.montage.M) || ...
            this.montage.Mind < 0
        % check montage index, if not good -> 0
        this.montage.Mind = 0;
    end
end

%-Check field order 
%--------------------------------------------------------------------------
fieldnames_order = {
    'type'
    'Nsamples'
    'Fsample'
    'timeOnset'
    'trials'
    'channels'
    'data'
    'fname'
    'path'
    'sensors'
    'fiducials'
    'transform'
    'condlist'
    'montage'
    'history'
    'other'};

[sel1, sel2] = match_str(fieldnames_order, fieldnames(this));
tempcell = struct2cell(this);
this = cell2struct(tempcell(sel2), fieldnames_order, 1);
