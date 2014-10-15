function D = spm_eeg_convert(S)
% Convert various M/EEG formats to SPM12 format
% FORMAT D = spm_eeg_convert(S)
% S                - string (filename) or struct (see below)
%
% If S is a struct it can have the optional following fields:
% S.dataset        - file name
% S.mode           - 'header'     - only convert the header without reading data
%                    'continuous' - convert data as continuous
%                    'epoched'    - convert data as epoched (requires data that is
%                                   already epoched or a trial definition in S.trl).
% S.timewin        - for continuous mode [start end] of data segment in sec (all if empty)
%                  - for epoched mode time window in PST ms
% S.outfile        - output file name (default 'spmeeg_' + input)
% S.channels       - 'all' - convert all channels
%                    or cell array of labels
% For epoched mode:
%
% S.trl            - [N x 3] trl matrix or name of the trial definition file
%                    containing 'trl' variable with such a matrix
% S.conditionlabels- labels for the trials in the data [default: 'Undefined']
%
%   or
%
% S.trialdef       - structure array for trial definition with fields
%     S.trialdef.conditionlabel - string label for the condition
%     S.trialdef.eventtype      - string
%     S.trialdef.eventvalue     - string, numeric or empty
%
%
% S.inputformat    - data type (optional) to force the use of specific data
%                    reader
% S.eventpadding   - the additional time period around each trial for which
%                    the events are saved with the trial (to let the user
%                    keep and use for analysis events which are outside
%                    trial borders), in seconds. [default: 0]
% S.blocksize      - size of blocks used internally to split large files
%                    [default: ~100Mb]
% S.checkboundary  - 1 - check if there are breaks in the file and do not
%                        read across those breaks [default]
%                    0 - ignore breaks (not recommended).
% S.saveorigheader - 1 - save original data header with the dataset
%                    0 - do not keep the original header [default]
%
% % D              - MEEG object (also written on disk)
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_convert.m 6244 2014-10-15 11:15:09Z vladimir $

SVNrev = '$Rev: 6244 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG convert'); spm('Pointer', 'Watch');

if ischar(S)
    temp      = S;
    S         = [];
    S.dataset = temp;
end

if ~isfield(S, 'dataset')
    error('Dataset must be specified.');
end

if ~isfield(S, 'outfile'),         S.outfile = ['spmeeg_' spm_file(S.dataset,'basename')]; end
if ~isfield(S, 'channels'),        S.channels = 'all';                                     end
if ~isfield(S, 'timewin'),         S.timewin   = [];                                       end
if ~isfield(S, 'blocksize'),       S.blocksize = 3276800;                                  end %100 Mb
if ~isfield(S, 'checkboundary'),   S.checkboundary = 1;                                    end
if ~isfield(S, 'eventpadding'),    S.eventpadding = 0;                                     end
if ~isfield(S, 'saveorigheader'),  S.saveorigheader = 0;                                   end
if ~isfield(S, 'conditionlabels'), S.conditionlabels = 'Undefined' ;                       end
if ~isfield(S, 'inputformat'),     S.inputformat = [] ;                                    end

if ~iscell(S.conditionlabels)
    S.conditionlabels = {S.conditionlabels};
end


if ~isfield(S, 'mode') || ~isequal(S.mode, 'header')
    % The header is read here in a recursive call and the above if avoids reading
    % it twice which might be expensive for some formats
    S1 = [];
    S1.mode           = 'header';
    S1.dataset        = S.dataset;
    S1.outfile        = S.outfile;
    S1.inputformat    = S.inputformat;
    Dhdr              = spm_eeg_convert(S1);
    hdr               = Dhdr.hdr;
    event             = Dhdr.events;
    eventsamples      = Dhdr.events(':', 'samples');
else
    %--------- Read and check header
    hdr = ft_read_header(S.dataset, 'headerformat', S.inputformat);
    
    if isfield(hdr, 'label')
        [unique_label,junk,ind] = unique(hdr.label);
        if length(unique_label) ~= length(hdr.label)
            warning(['Data file contains several channels with ',...
                'the same name. These channels cannot be processed and will be disregarded']);
            % This finds the repeating labels and removes all their occurences
            sortind     = sort(ind);
            [junk,ind2] = setdiff(hdr.label, unique_label(sortind(diff(sortind)==0)));
            hdr.label   = hdr.label(ind2);
            hdr.nChans  = length(hdr.label);
        end
    end
    
    
    %--------- Read and prepare events
    
    try
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
                bytes  = bytes(:, end:-1:end-7);
                val    = num2cell(bin2dec(bytes));
                [event(ind).value] = deal(val{:});
            end
        end
        
    catch
        warning(['Could not read events from file ' S.dataset]);
        event = [];
    end
    
    % Replace samples with time
    if numel(event)>0
        for i = 1:numel(event)
            event(i).time = event(i).sample./hdr.Fs;
        end
    end
end


if ~isfield(S, 'mode')
    if (hdr.nTrials == 1)
        S.mode = 'continuous';
    else
        S.mode = 'epoched';
    end
end

%--------- Start making the header

D = [];
D.Fsample = hdr.Fs;

%--------- Select channels

if ~strcmp(S.channels, 'all') %Dhdr should be available in this case
    chansel = selectchannels(Dhdr, S.channels);
else
    if isfield(hdr, 'nChans')
        chansel = 1:hdr.nChans;
    else
        chansel = 1:length(hdr.label);
    end
end

nchan = length(chansel);

D.channels = repmat(struct('bad', 0), 1, nchan);

if isfield(hdr, 'label')
    [D.channels(:).label] = deal(hdr.label{chansel});
end
%--------- Preparations specific to reading mode (continuous/epoched)

if ismember(S.mode, {'continuous', 'header'})
    
    if isempty(S.timewin)
        if hdr.nTrials == 1
            segmentbounds = [1 hdr.nSamples];
        elseif ~S.checkboundary || isequal(S.mode, 'header')
            segmentbounds = [1 hdr.nSamples*hdr.nTrials];
        else
            error('The data cannot be read without ignoring trial borders');
        end
        timewindow = segmentbounds./D.Fsample;
    else
        timewindow = S.timewin;
        segmentbounds = round(timewindow.*D.Fsample);
        segmentbounds(1) = max(segmentbounds(1), 1);
    end
    
    
    %--------- Sort events and put in the trial
    
    if ~isempty(event)
        try
            event = rmfield(event, {'offset', 'sample'});
        end
        event = select_events(event, ...
            [timewindow(1)-S.eventpadding timewindow(2)+S.eventpadding]);
    end
    
    D.trials.label  = S.conditionlabels{1};
    D.trials.events = event;
    D.trials.onset  = timewindow(1);
    
    %--------- Break too long segments into blocks
    
    nblocksamples = floor(S.blocksize/nchan);
    nsampl = diff(segmentbounds)+1;
    
    trl = segmentbounds(1):nblocksamples:segmentbounds(2);
    if (trl(end)==segmentbounds(2))
        trl = trl(1:(end-1));
    end
    
    trl = [trl(:) [trl(2:end)-1 segmentbounds(2)]'];
    
    ntrial = size(trl, 1);
    
    readbytrials = 0;
    
    
    D.timeOnset = (trl(1,1)-1)./hdr.Fs;
 
    D.Nsamples = nsampl;
else % Read by trials
    if isfield(S, 'trl') || isfield(S, 'trialdef')
        if isfield(S, 'trl')
            if ischar(S.trl)
                trl = getfield(load(S.trl, 'trl'), 'trl');
                conditionlabels = getfield(load(S.trl, 'conditionlabels'), 'conditionlabels');
            else
                trl = S.trl;
                conditionlabels = S.conditionlabels;
            end
        else            
            S1          = [];
            S1.D        = Dhdr;
            S1.timewin  = S.timewin;
            S1.trialdef = S.trialdef;
            S1.reviewtrials = 0;
            S1.save  = 0;
            [trl, conditionlabels] = spm_eeg_definetrial(S1);
        end
        
        trl = double(trl);
        
        if size(trl, 2) >= 3
            D.timeOnset = unique(trl(:, 3))./D.Fsample;
            trl = trl(:, 1:2);
        else
            D.timeOnset = 0;
        end
        
        if length(D.timeOnset) > 1
            error('All trials should have identical baseline');
        end
        
        if ~iscell(conditionlabels)
            conditionlabels = {conditionlabels};
        end
        
        if numel(conditionlabels) == 1
            conditionlabels = repmat(conditionlabels, 1, size(trl, 1));
        end
        
        readbytrials = 0;
    else
        try
            trialind = sort([strmatch('trial', {event.type}, 'exact'), ...
                strmatch('average', {event.type}, 'exact')]);
            trl = [eventsamples(trialind).sample];
            trl = double(trl(:));
            trl = [trl  trl+double([event(trialind).duration]')-1];
            
            try
                offset = unique([event(trialind).offset]);
            catch
                offset = [];
            end
            
            if length(offset) == 1 && offset~=0
                D.timeOnset = offset/D.Fsample;
            else            
                D.timeOnset = -hdr.nSamplesPre/hdr.Fs;
            end
            conditionlabels = {};
            for i = 1:length(trialind)
                if isempty(event(trialind(i)).value)
                    conditionlabels{i} = S.conditionlabels{1};
                else
                    if all(ischar(event(trialind(i)).value))
                        conditionlabels{i} = event(trialind(i)).value;
                    else
                        conditionlabels{i} = num2str(event(trialind(i)).value);
                    end
                end
            end
            if  hdr.nTrials>1 && size(trl, 1)~=hdr.nTrials
                warning('Mismatch between trial definition in events and in data. Ignoring events');
                readbytrials = 1;
            else
                readbytrials = 0;
            end
            
            event = event(setdiff(1:numel(event), trialind));
        catch
            if hdr.nTrials == 1
                error('Could not define trials based on data. Use continuous option or trial definition file.');
            else
                readbytrials = 1;
            end
        end
    end
    if readbytrials
        nsampl = hdr.nSamples;
        ntrial = hdr.nTrials;
        trl = zeros(ntrial, 2);
        if exist('conditionlabels', 'var') ~= 1 || length(conditionlabels) ~= ntrial
            conditionlabels = repmat(S.conditionlabels, 1, ntrial);
        end
    else
        nsampl = unique(diff(trl, [], 2))+1;
        if length(nsampl) > 1
            error('All trials should have identical lengths');
        end
        
        inbounds = (trl(:,1)>=1 & trl(:, 2)<=hdr.nSamples*hdr.nTrials)';
        
        rejected = find(~inbounds);
        
        if ~isempty(rejected)
            trl = trl(inbounds, :);
            conditionlabels = conditionlabels(inbounds);
            warning([S.dataset ': Trials ' num2str(rejected) ' not read - out of bounds']);
        end
        
        ntrial = size(trl, 1);
        
        if ntrial == 0
            warning([S.dataset ': No trials to read. Bailing out.']);
            D = [];
            return;
        end
    end
    D.Nsamples = nsampl;
    if isfield(event, 'sample')
        event = rmfield(event, 'sample');
    end
end

%--------- Prepare for reading the data
outpath = spm_file(S.outfile,'fpath');
outfile = spm_file(S.outfile,'basename');
if isempty(outfile), outfile = 'spm8'; end

D.path = outpath;
D.fname = [outfile '.mat'];

if ~isequal(S.mode, 'header')
    
    % These are the units used for channel data in SPM
    chanunit  = units(Dhdr, chansel);
    chanunit(strcmp('T',   chanunit)) = {'fT'};
    chanunit(strcmp('T/m', chanunit)) = {'fT/mm'};
    chanunit(strcmp('V',   chanunit)) = {'uV'};
    
    if isequal(S.mode, 'continuous')
        D.data = file_array(fullfile(D.path, [outfile '.dat']), [nchan nsampl], 'float32-le');
    else
        D.data = file_array(fullfile(D.path, [outfile '.dat']), [nchan nsampl ntrial], 'float32-le');
    end
    
    % physically initialise file
    initialise(D.data);
    
    spm_progress_bar('Init', ntrial, 'reading and converting'); drawnow;
    if ntrial > 100, Ibar = floor(linspace(1, ntrial,100));
    else Ibar = 1:ntrial; end
    
    %--------- Read the data
    
    offset = 1;
    for i = 1:ntrial
        spm_progress_bar('Set','ylabel','reading...');
        if readbytrials
            dat = ft_read_data(S.dataset,'header',  hdr, 'begtrial', i, 'endtrial', i,...
                'chanindx', chansel, 'chanunit', chanunit, 'checkboundary', S.checkboundary, 'dataformat', S.inputformat);
        else
            dat = ft_read_data(S.dataset,'header',  hdr, 'begsample', trl(i, 1), 'endsample', trl(i, 2),...
                'chanindx', chansel, 'chanunit', chanunit, 'checkboundary', S.checkboundary, 'dataformat', S.inputformat); 
        end
        
        % Sometimes ft_read_data returns sparse output
        dat = full(dat);
        
        spm_progress_bar('Set','ylabel','writing...');
        switch S.mode
            case 'continuous'
                nblocksamples = size(dat,2);
                
                D.data(:, offset:(offset+nblocksamples-1)) = dat;
                
                offset = offset+nblocksamples;
            case 'epoched'
                D.data(:, :, i) = dat;
                D.trials(i).label = conditionlabels{i};
                D.trials(i).onset = trl(i, 1)./D.Fsample;
                D.trials(i).events = select_events(event, ...
                    [ trl(i, 1)./D.Fsample-S.eventpadding  trl(i, 2)./D.Fsample+S.eventpadding]);
        end
        
        if ismember(i, Ibar)
            spm_progress_bar('Set', i);
        end
        
    end
    
    spm_progress_bar('Clear');
    
else
    D.data = [];
end

% Specify sensor positions and fiducials
if isfield(hdr, 'grad')
    D.sensors.meg = hdr.grad;
end

if isfield(hdr, 'elec')
    elec = hdr.elec;
else
    try
        elec = ft_read_sens(S.dataset, 'fileformat', S.inputformat);
        % It might be that read_sens will return the grad for MEG datasets
        if any(ismember({'ori', 'coilori', 'coilpos'}, fieldnames(elec)))
            elec = [];
        end
    catch
        elec = [];
    end
end

if ~isempty(elec)   
    D.sensors.eeg = elec;
else
    sw = warning('off','backtrace');
    warning('Could not obtain electrode locations automatically.');
    warning(sw);
end


try
    D.fiducials = ft_convert_units(ft_read_headshape(S.dataset, 'fileformat', S.inputformat), 'mm');
catch
    sw = warning('off','backtrace');
    warning('Could not obtain fiducials automatically.');
    warning(sw);
end

%--------- Create meeg object
D = meeg(D);

% history
D = D.history('spm_eeg_convert', S);

if isfield(hdr, 'orig')
    if S.saveorigheader
        D.origheader = hdr.orig;
    end
    
    % Uses fileio function to get the information about channel types stored in
    % the original header.
    origchantypes = ft_chantype(hdr);
    [sel1, sel2] = spm_match_str(D.chanlabels, hdr.label);
    origchantypes = origchantypes(sel2);
    if length(strmatch('unknown', origchantypes, 'exact')) ~= numel(origchantypes)
        D.origchantypes = struct([]);
        D.origchantypes(1).label = hdr.label(sel2);
        D.origchantypes(1).type = origchantypes;
    end
end

S1 = [];
S1.task = 'defaulttype';
S1.D = D;
S1.updatehistory = 0;
D = spm_eeg_prep(S1);

% Assign default EEG sensor positions if possible
if ~isempty(strmatch('EEG', D.chantype, 'exact'))
    if isempty(D.sensors('EEG'))
        S1 = [];
        S1.task = 'defaulteegsens';
        S1.updatehistory = 0;
        S1.D = D;
        
        D = spm_eeg_prep(S1);
    else
        S1 = [];
        S1.task = 'project3D';
        S1.modality = 'EEG';
        S1.updatehistory = 0;
        S1.D = D;
        
        D = spm_eeg_prep(S1);
    end
end

% Create 2D positions for MEG
% by projecting the 3D positions to 2D
if ~isempty(strmatch('MEG', D.chantype)) && ~isempty(D.sensors('MEG'))
    S1 = [];
    S1.task = 'project3D';
    S1.modality = 'MEG';
    S1.updatehistory = 0;
    S1.D = D;
    
    D = spm_eeg_prep(S1);
end

% If channel units are available, store them.
if isequal(S.mode, 'header') 
    [sel1, sel2] = spm_match_str(D.chanlabels, hdr.label);
    D = units(D, sel1, hdr.chanunit(sel2));
else
    D = units(D, ':', chanunit);
end

% The conditions will later be sorted in the original order they were defined.
if isfield(S, 'trialdef')
    D = condlist(D, {S.trialdef(:).conditionlabel});
end

% This is for the recursive call to work properly
if isequal(S.mode, 'header')
    D.hdr = hdr;
end

if ~isequal(S.mode, 'header')
    save(D);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG convert: done'); spm('Pointer', 'Arrow');


%==========================================================================
% select_events
%==========================================================================
function event = select_events(event, timeseg)
% Utility function to select events according to time segment
% FORMAT event = select_events(event, timeseg)
if iscell(event)
    event = event{1};
end

if ~isempty(event)
    [time,ind] = sort([event(:).time]);
    
    selectind  = ind(time>=timeseg(1) & time<=timeseg(2));
    
    event      = event(selectind);
end
