function spm_eeg_convert(S)
% Main function for converting different M/EEG formats to SPM8 format.
% FORMAT spm_eeg_convert(S)
% S - can be string (file name) or struct.
%
% If S is a struct it can have the following fields:
% S.dataset - file name
% S.continuous - 1 - convert data as continuous
%                0 - convert data as epoched (requires data that is already
%                    epoched or a trial definition file).
% S.timewindow - [start end] in sec. Boundaries for a sub-segment of
%                continuous data (default - all).
% S.outfile - name base for the output files (default - the same as input)
% S.channels - 'all' - convert all channels
%               cell array of labels
% S.chanfile - name of the channel selection file
% S.usetrials - 1 - take the trials as defined in the data (default)
%               0 - use trial definition file even though the data is
%                   already epoched.
% S.trlfile - name of the trial definition file
% S.datatype - data type for the data file one of
%              'int16','int32','float32' (default), 'float64'
% S.eventpadding - in sec - the additional time period around each trial
%               for which the events are saved with the trial (to let the
%               user keep and use for analysis events which are outside
%               trial borders). Default - 0.
% S.conditionlabels - labels for the trials in the data Default - 'Undefined'
% S.blocksize - size of blocks used internally to split large files
%               default ~100Mb.
% S.checkboundary - 1 - check if there are breaks in the file and do not read
%                       across those breaks (default).
%                   0 - ignore breaks (not recommended).
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_convert.m 1548 2008-05-06 15:07:31Z vladimir $

[Finter] = spm('FnUIsetup','MEEG data conversion ',0);

if ischar(S)
    temp = S;
    S=[];
    S.dataset = temp;
end

if ~isfield(S, 'dataset')
    error('Dataset must be specified.');
end

if ~isfield(S, 'outfile'),         S.outfile = ['spm8_' spm_str_manip(S.dataset,'tr')];     end
if ~isfield(S, 'channels'),        S.channels = 'all';                                      end
if ~isfield(S, 'timewindow'),      S.timewindow = [];                                       end
if ~isfield(S, 'blocksize'),       S.blocksize = 3276800;                                   end  %100 Mb
if ~isfield(S, 'checkboundary'),   S.checkboundary = 1;                                     end
if ~isfield(S, 'usetrials'),       S.usetrials = 1;                                         end
if ~isfield(S, 'datatype'),        S.datatype = 'float32';                                  end
if ~isfield(S, 'eventpadding'),    S.eventpadding = 0;                                      end
if ~isfield(S, 'conditionlabel'),  S.conditionlabel = 'Undefined';                          end

%--------- Read and check header

hdr = fileio_read_header(S.dataset, 'fallback', 'biosig');

if isfield(hdr, 'label')
    [unique_label junk ind]=unique(hdr.label);
    if length(unique_label)~=length(hdr.label)
        warning(['Data file contains several channels with ',...
            'the same name. These channels cannot be processed and will be disregarded']);
        % This finds the repeating labels and removes all their occurences
        sortind=sort(ind);
        [junk ind2]=setdiff(hdr.label, unique_label(sortind(find(diff(sortind)==0))));
        hdr.label=hdr.label(ind2);
        hdr.nChans=length(hdr.label);
    end
end

if ~isfield(S, 'continuous')
    S.continuous = (hdr.nTrials == 1);
end

%--------- Read and prepare events

try
    event = fileio_read_event(S.dataset, 'detectflank', 'both');
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

%--------- Start making the header

D = [];
D.Fsample = hdr.Fs;

%--------- Select channels

if ~strcmp(S.channels, 'all')
    [junk, chansel] = spm_match_str(S.channels, hdr.label);
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

if S.continuous

    if isempty(S.timewindow)
        segmentbounds = [1 hdr.nSamples];
        S.timewindow = segmentbounds./D.Fsample;
    else
        segmentbounds = S.timewindow.*D.Fsample;
        segmentbounds(1) = max(segmentbounds(1), 1);
    end

    %--------- Sort events and put in the trial

    if ~isempty(event)
        event = rmfield(event, {'offset', 'sample'});
        event = select_events(event, ...
            [S.timewindow(1)-S.eventpadding S.timewindow(2)+S.eventpadding]);
    end

    D.trials.label = S.conditionlabel;
    D.trials.events = event;
    D.trials.onset = S.timewindow(1);

    %--------- Break too long segments into blocks

    nblocksamples = floor(S.blocksize/nchan);
    nsampl = diff(segmentbounds)+1;

    trl = [segmentbounds(1):nblocksamples:segmentbounds(2)];
    if (trl(end)==segmentbounds(2))
        trl = trl(1:(end-1));
    end

    trl = [trl(:) [trl(2:end)-1 segmentbounds(2)]'];

    ntrial = size(trl, 1);

    readbytrials = 0;

    D.timeOnset = (-hdr.nSamplesPre+trl(1,1)-1)./hdr.Fs;
    D.Nsamples = nsampl;
else % Read by trials
    if ~S.usetrials
        if ~isfield(S, 'trl')
            trl = getfield(load(S.trlfile, 'trl'), 'trl');
        else
            trl = S.trl;
        end
        if size(trl, 2) >= 3
            D.timeOnset = unique(trl(:, 3))./D.Fsample;
            trl = trl(:, 1:2);
        else
            D.timeOnset = 0;
        end

        if length(D.timeOnset) > 1
            error('All trials should have identical baseline');
        end
        if isfield(S, 'conditionlabels')
            conditionlabels = S.conditionlabels;
        else
            conditionlabels = getfield(load(S.trlfile, 'conditionlabels'), 'conditionlabels');
        end

        if numel(conditionlabels) == 1
            conditionlabels = repmat(conditionlabels, 1, size(trl, 1));
        end

        readbytrials = 0;
    else
        try
            trialind = find(strcmp('trial', {event.type}));
            trl = [event(trialind).sample];
            trl = trl(:);
            trl = [trl  trl+[event(trialind).duration]'-1];

            try
                offset = unique([event(trialind).offset]);
            catch
                offset = [];
            end
            if length(offset) == 1
                D.timeOnset = offset/D.Fsample;
            else
                D.timeOnset = 0;
            end
            conditionlabels = {};
            for i = 1:length(trialind)
                if isempty(event(trialind(i)).value)
                    conditionlabels{i} = S.conditionlabel;
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
            end

            readbytrials = 1;
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
        if length(conditionlabels) ~= ntrial
            conditionlabels = repmat({S.conditionlabel}, 1, ntrial);
        end
    else
        nsampl = unique(diff(trl, [], 2))+1;
        if length(nsampl) > 1
            error('All trials should have identical lengths');
        end

        inbounds = (trl(:,1)>1 & trl(:, 2)<=hdr.nSamples)';

        rejected = find(~inbounds);

        if ~isempty(rejected)
            trl = trl(find(inbounds), :);
            warning([S.dataset ': Trials ' num2str(rejected) ' not read - out of bounds']);
        end

        ntrial = size(trl, 1);
    end
    D.Nsamples = nsampl;
    if isfield(event, 'sample')
        event = rmfield(event, 'sample');
    end
end

%--------- Prepare for reading the data
D.data.fnamedat = [S.outfile '.dat'];
D.data.datatype = S.datatype;

if S.continuous
    datafile = file_array(D.data.fnamedat, [nchan nsampl], S.datatype);
else
    datafile = file_array(D.data.fnamedat, [nchan nsampl ntrial], S.datatype);
end

% physically initialise file
datafile(end,end) = 0;

spm('Pointer', 'Watch');drawnow;

spm_progress_bar('Init', ntrial, 'reading and converting'); drawnow;
if ntrial > 100, Ibar = floor(linspace(1, ntrial,100));
else Ibar = [1:ntrial]; end

%--------- Read the data

offset = 1;
for i = 1:ntrial
    if readbytrials
        dat = fileio_read_data(S.dataset,'header',  hdr, 'begtrial', i, 'endtrial', i,...
            'chanindx', chansel, 'checkboundary', S.checkboundary, 'fallback', 'biosig');
    else
        dat = fileio_read_data(S.dataset,'header',  hdr, 'begsample', trl(i, 1), 'endsample', trl(i, 2),...
            'chanindx', chansel, 'checkboundary', S.checkboundary, 'fallback', 'biosig');
    end

    % Sometimes fileio_read_data returns sparse output
    dat = full(dat);

    if S.continuous
        nblocksamples = size(dat,2);

        datafile(:, offset:(offset+nblocksamples-1)) = dat;

        offset = offset+nblocksamples;
    else
        datafile(:, :, i) = dat;
        D.trials(i).label = conditionlabels{i};
        D.trials(i).onset = trl(i, 1)./D.Fsample;
        D.trials(i).events = select_events(event, ...
            [ trl(i, 1)./D.Fsample-S.eventpadding  trl(i, 2)./D.Fsample+S.eventpadding]);
    end

    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end

end

spm_progress_bar('Clear');

spm('Pointer', 'Arrow');drawnow;

% Specify sensor positions and fiducials
if isfield(hdr, 'grad')
    D.sensors.meg = forwinv_convert_units(hdr.grad, 'mm');
elseif isfield(hdr, 'elec')
    D.sensors.eeg = forwinv_convert_units(hdr.elec, 'mm');
else
    try
        D.sensors.eeg = forwinv_convert_units(fileio_read_sens(S.dataset), 'mm');
    catch
        warning('Could not obtain electrode locations automatically.');
    end
end

try
    D.fiducials = forwinv_convert_units(fileio_read_headshape(S.dataset), 'mm');
catch
    warning('Could not obtain fiducials automatically.');
end

%--------- Create meeg object
D.history(1).fun = 'spm_eeg_convert';
D.history(1).args = {S};
D.fname = [S.outfile '.mat'];

D = meeg(D);

% Set channel types to default
D = chantype(D, [], []);

if S.continuous
    D = type(D, 'continuous');
else
    D = type(D, 'single');
end

save(D);

function event = select_events(event, timeseg)
% Utility function to select events according to time segment
% FORMAT event = select_events(event, timeseg)

if ~isempty(event)
    [time ind] = sort([event(:).time]);

    selectind = ind(time>=timeseg(1) & time<=timeseg(2));

    event = event(selectind);
end