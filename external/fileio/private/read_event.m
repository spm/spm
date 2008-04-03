function [event] = read_event(filename, varargin)

% READ_EVENT reads all events from an EEG/MEG dataset and returns
% them in a well defined structure. It is a wrapper around different
% EEG/MEG file importers, directly supported formats are CTF, Neuromag,
% EEP, BrainVision, Neuroscan and Neuralynx.
%
% Use as
%   [event] = read_event(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'eventformat'   string
%   'header'        structure, see READ_HEADER
% Furthermore, you can specify optional arguments as key-value pairs for
% filtering the events, e.g. to select only events of a specific type. See
% FILTER_EVENT for more details.
%
% This function returns an event structure with the following fields
%   event.type      string
%   event.sample    expressed in samples, first sample of file is 1
%   event.value     number or string
%   event.offset    expressed in samples
%   event.duration  expressed in samples
%
% The event type and sample fields are always defined, other fields can be empty,
% depending on the type of event file. Events are sorted by the sample on
% which they occur. After reading the event structure, you can use the
% following tricks to extract information about those events in which you
% are interested.
%
% Determine the different event types
%   unique({event.type})
%
% Get the index of all trial events
%   find(strcmp('trial', {event.type}))
%
% Make a vector with all triggers that occurred on the backpanel
%   [event(find(strcmp('backpanel trigger', {event.type}))).value]
%
% Find the events that occurred in trial 26
%   t=26; samples_trials = [event(find(strcmp('trial', {event.type}))).sample];
%   find([event.sample]>samples_trials(t) & [event.sample]<samples_trials(t+1))
%
% See also READ_HEADER, READ_DATA, WRITE_DATA, WRITE_EVENT, FILTER_EVENT

% Copyright (C) 2004-2007, Robert Oostenveld
%
% $Log: read_event.m,v $
% Revision 1.51  2007/12/20 19:06:57  roboos
% Added filtering base on event number (minnumber and maxnumber), implemented in low level for neuralynx_nev and for the rest of the formats in filter event. If event.number is not present everything still should work as it used to.
%
% Revision 1.50  2007/12/19 15:24:37  roboos
% fixed typo, ttl should be bin
%
% Revision 1.49  2007/12/19 11:25:59  roboos
% fixed typo, added "(", thanks to Mahdi
%
% Revision 1.48  2007/12/19 09:29:29  roboos
% implmenented events for neuralynx_sdma file, and merged it with the dma, ttl and bin formats
% cleaned up the consistent handling of the ttl values for dma, ttl and bin formats
%
% Revision 1.47  2007/12/18 16:57:35  roboos
% use value=trigger instead of ttl (which was a change in the previous commit) to avoid breaking Thilo's scripts
%
% Revision 1.46  2007/12/18 16:51:34  roboos
% added some filtering options
% some minor changes to neuralynx_nev, related to the changes in teh low level function
% added support for neuralynx_bin, which is treated similar as ttl and dma
%
% Revision 1.45  2007/12/17 12:59:52  roboos
% reimplemented the event detection for ttl and dma after discussion with Thilo, also merged the two implementations
%
% Revision 1.44  2007/12/17 08:24:28  roboos
% added support for nexstim_nxe, thanks to Vladimir
% the low-level code has not been tested by myself
%
% Revision 1.43  2007/12/12 16:51:07  roboos
% made a faster implementation for a nev file inside a neuralynx_ds dataset directory, but strings are not supported any more
%
% Revision 1.42  2007/12/12 11:29:06  roboos
% chedk for presence of timestamp prior to trying to concatenate events
%
% Revision 1.41  2007/12/12 11:10:55  roboos
% moved declaration of global variable to the begin of the function
%
% Revision 1.40  2007/12/12 11:09:22  roboos
% added selective reading of events for some files (dma, ttl, bdf, 4d, ctf partially)
%
% Revision 1.39  2007/11/07 10:49:06  roboos
% cleaned up the reading and writing from/to mysql database, using db_xxx helper functions (see mysql directory)
%
% Revision 1.38  2007/11/05 17:02:07  roboos
% some cosmetic changes, nothing functional
%
% Revision 1.37  2007/10/16 12:34:12  roboos
% use recursion to read from multiple event sources
% implemented fcdc_global
%
% Revision 1.36  2007/10/02 16:07:13  roboos
% fixed conversion of negative status values and reasignment of bit24 (thanks to Philip)
%
% Revision 1.35  2007/10/02 09:28:06  roboos
% removed teh double rtepresentatino of trigger for bdf, since STATUS is also good enough
%
% Revision 1.34  2007/10/02 09:13:43  roboos
% biosemi_bdf: make sure that the sign bit is propperly re-inserted as bit24 for the status channel
%
% Revision 1.33  2007/10/01 13:43:21  roboos
% reimplemented the biosemi bdf trigger detection, now also for the status bits
%
% Revision 1.32  2007/09/13 09:49:34  roboos
% moved declaration of persistent variale to beginning
% added inactive piece of code for bdf (see NICI version)
%
% Revision 1.31  2007/08/21 17:00:57  chrhes
% updated some documentation, removed the commented-out section of code to do
% with the event filtering options; these are used "as is" in the call to
% FILTER_EVENT at the end of the function
%
% Revision 1.30  2007/08/01 12:24:22  roboos
% added filename as argument to read_shm_event, added comments, disabled unused 
% keyval filtering arguments
%
% Revision 1.29  2007/08/01 09:57:01  roboos
% moved all code related to ctf shared memory to seperate functions
%
% Revision 1.28  2007/07/30 12:14:59  roboos
% updated documentation for filtering
% implemented filtering for ctf_shm
%
% Revision 1.27  2007/07/27 12:17:20  roboos
% fixed big in ctf raw trigger channel, which caused a trigger with value=1 not 
% to be detected implemented support for ctf_shm
% reuse the file header if specified as optional input argument and do not read 
% again
%
% Revision 1.26  2007/07/04 13:20:51  roboos
% added support for egi_egis/egia, thanks to Joseph Dien
%
% Revision 1.25  2007/06/13 13:33:54  roboos
% changed the mysql code to reflect the updated event table structure
% removed type and subtype from the insert query
% added a call to filter_event at the end of the function
%
% Revision 1.24  2007/06/13 09:57:32  roboos
% only test for presence of fields if event is not empty
% fixed bug in mysql, in case event table is empty
%
% Revision 1.23  2007/06/13 08:06:21  roboos
% updated help
%
% Revision 1.22  2007/06/12 19:35:52  roboos
% implemented support for reading events from mysql database
%
% Revision 1.21  2007/06/11 13:52:24  roboos
% split the event reading from neuralynx_dma and neuralynx_ttl
% read timestamps from tsl/tsh files, optionally use pre-specified header to 
% correct the sample numbers
%
% Revision 1.20  2007/06/07 12:44:28  chrhes
% updated some documentation
%
% Revision 1.19  2007/06/06 21:55:37  chrhes
% fixed a small bug to do with a string comparison
%
% Revision 1.18  2007/06/06 18:19:07  chrhes
% added initial implementation of reading events from the serial port
%
% Revision 1.17  2007/06/06 07:12:48  roboos
% switched to using filetype_check_uri for detection and parsing of filename
%
% Revision 1.16  2007/05/31 09:53:21  roboos
% implemented reading events from a plain matlab file
%
% Revision 1.15  2007/05/31 09:15:13  roboos
% added placeholder for tcpsocket
%
% Revision 1.14  2007/05/15 15:01:44  roboos
% changed handling of the seperate brainvision header and marker file
%
% Revision 1.13  2006/12/13 15:40:10  roboos
% renamed Parallel_in into ttl for neuralynx_dma
% renamed the function read_neuralynx_event into read_neuralynx_nev (consistent 
% with the file extension)
%
% Revision 1.12  2006/12/04 10:37:27  roboos
% added support for ns_avg
%
% Revision 1.11  2006/09/18 21:51:38  roboos
% implemented support for fcdc_matbin, i.e. a dataset consisting of a matlab 
% file with header and events and a seperate binary datafile
%
% Revision 1.10  2006/09/18 14:22:54  roboos
% implemented support for 4D-BTi dataformat
%
% Revision 1.9  2006/08/28 10:13:03  roboos
% use seperate filetype_check_extension function instead of subfunction, removed 
% subfunction
%
% Revision 1.8  2006/07/26 07:51:33  roboos
% fixed bug for brainvision_vmrk in case second field is empty (thanks to 
% Stephan Bickel)
%
% Revision 1.7  2006/06/26 08:46:37  roboos
% read stim channels from CTF as continuous
%
% Revision 1.6  2006/06/22 07:54:34  roboos
% do not read the complete data for ns_cnt but only the header (includes the 
% events), thanks to Gijs
%
% Revision 1.5  2006/06/20 11:23:00  roboos
% fixed bug for ctf, sensSype was moved to hdr.orig
%
% Revision 1.4  2006/06/19 10:32:16  roboos
% added documentation
%
% Revision 1.3  2006/06/19 08:14:17  roboos
% updated documentation
%
% Revision 1.2  2006/06/07 10:16:47  roboos
% changed one occurence of read_fcdc_header into read_header
%
% Revision 1.1  2006/06/07 09:32:20  roboos
% new implementation based on the read_fcdc_xxx functions, now with
% variable (key-val) input arguments, changed the control structure
% in the rpobram (switch instead of ifs), allow the user to specify
% the file format, allow the user to specify either a sample selection
% or a block selection. The reading functionality should not have
% changed compared to the read_fcdc_xxx versions.
%

global event_queue        % for fcdc_global
persistent db_blob        % for fcdc_mysql
if isempty(db_blob)
  db_blob = 0;
end

if iscell(filename)
  % use recursion to read from multiple event sources
  event = [];
  for i=1:numel(filename)
    tmp   = read_event(filename{i}, varargin{:});
    event = catevent(event(:), tmp(:));
  end
  return
end
    
% get the options
eventformat      = keyval('eventformat', varargin);
hdr              = keyval('header',      varargin);
% this allows to read only events in a certain range, supported for selected data formats only
flt_type         = keyval('type',         varargin);  
flt_value        = keyval('value',        varargin);  
flt_minsample    = keyval('minsample',    varargin);  
flt_maxsample    = keyval('maxsample',    varargin);
flt_mintimestamp = keyval('mintimestamp', varargin);  
flt_maxtimestamp = keyval('maxtimestamp', varargin);
flt_minnumber    = keyval('minnumber', varargin);  
flt_maxnumber    = keyval('maxnumber', varargin);

% determine the filetype
if isempty(eventformat)
  eventformat = filetype(filename);
end

switch eventformat
  case 'brainvision_vhdr'
    % read the headerfile belonging to the dataset and try to determine the corresponding markerfile
    eventformat = 'brainvision_vmrk';
    hdr = read_brainvision_vhdr(filename);
    % replace the filename with the filename of the markerfile
    if ~isfield(hdr, 'MarkerFile') || isempty(hdr.MarkerFile)
      filename = [];
    else
      [p, f, e] = fileparts(filename);
      filename = fullfile(p, hdr.MarkerFile);
    end
end

% start with an empty event structure
event = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the events with the low-level reading function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch eventformat

  case 'fcdc_global'
    event = event_queue;

  case {'4d_pdf', '4d_m4d', '4d_xyz'}
    if isempty(hdr)
      hdr = read_header(filename);
    end
    % add the trials to the event structure
    for i=1:hdr.nTrials
      event(end+1).type     = 'trial';
      event(end  ).sample   = (i-1)*hdr.nSamples + 1;
      event(end  ).value    = [];
      event(end  ).offset   = -hdr.nSamplesPre;
      event(end  ).duration = hdr.nSamples;
    end
    % specify the range to search for triggers, default is the complete file
    if ~isempty(flt_minsample)
      begsample = flt_minsample;
    else
      begsample = 1;
    end
    if ~isempty(flt_maxsample)
      endsample = flt_maxsample;
    else
      endsample = hdr.nSamples*hdr.nTrials;
    end
    % read the data from the trigger channel
    dat = read_data(filename, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', hdr.orig.TriggerIndex);
    for i=1:size(dat,1)
      % detect the flanks in the trigger channel
      % FIXME the trigger channel looks a bit sloppy, i.e. not really sharp transitions
      flank = [0 diff(dat(i,:), 1, 2)~=0];
      indx = find(flank);
      trig = dat(indx);
      for j=1:length(trig)
        event(end+1).type     = hdr.label{hdr.orig.TriggerIndex(i)};
        event(end  ).sample   = indx(j) + begsample - 1;
        event(end  ).value    = trig(j);
        event(end  ).offset   = 0;
        event(end  ).duration = 0;
      end
    end

  case {'besa_avr', 'besa_swf'}
    if isempty(hdr)
      hdr = read_header(filename);
    end
    event(end+1).type     = 'average';
    event(end  ).sample   = 1;
    event(end  ).duration = hdr.nSamples;
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).value    = [];

  case {'biosemi_bdf', 'bham_bdf'}
    % read the header, required to determine the stimulus channels and trial specification
    if isempty(hdr)
      hdr = read_header(filename);
    end
    % specify the range to search for triggers, default is the complete file
    if ~isempty(flt_minsample)
      begsample = flt_minsample;
    else
      begsample = 1;
    end
    if ~isempty(flt_maxsample)
      endsample = flt_maxsample;
    else
      endsample = hdr.nSamples*hdr.nTrials;
    end
    % find the STATUS channel and read the values from it
    schan = find(strcmp(upper(hdr.label),'STATUS'));
    sdata = read_data(filename, hdr, begsample, endsample, schan);

    % find indices of negative numbers
    bit24i = find(sdata < 0);
    % make number positive and preserve bits 0-22
    sdata(bit24i) = bitcmp(abs(sdata(bit24i))-1,24);
    % re-insert the sign bit on its original location, i.e. bit24 
    sdata(bit24i) = sdata(bit24i)+(2^(24-1));
    % typecast the data to ensure that the status channel is represented in 32 bits
    sdata = uint32(sdata);        
    
    byte1 = 2^8  - 1;
    byte2 = 2^16 - 1 - byte1;
    byte3 = 2^24 - 1 - byte1 - byte2;

    % get the respective status and trigger bits
    trigger = bitand(sdata, bitor(byte1, byte2)); %  contained in the lower two bytes
    epoch   = int8(bitget(sdata, 16+1));
    cmrange = int8(bitget(sdata, 20+1));
    battery = int8(bitget(sdata, 22+1));

    % determine when the respective status bits go up or down
    flank_trigger = diff([0 trigger]);
    flank_epoch   = diff([0 epoch ]);
    flank_cmrange = diff([0 cmrange]);
    flank_battery = diff([0 battery]);

    for i=find(flank_trigger>0)
      event(end+1).type   = 'STATUS';
      event(end  ).sample = i + begsample - 1;
      event(end  ).value  = double(trigger(i));
    end

    for i=find(flank_epoch==1)
      event(end+1).type   = 'Epoch';
      event(end  ).sample = i;
    end

    for i=find(flank_cmrange==1)
      event(end+1).type   = 'CM_in_range';
      event(end  ).sample = i;
    end

    for i=find(flank_cmrange==-1)
      event(end+1).type   = 'CM_out_of_range';
      event(end  ).sample = i;
    end

    for i=find(flank_battery==1)
      event(end+1).type   = 'Battery_low';
      event(end  ).sample = i;
    end

    for i=find(flank_battery==-1)
      event(end+1).type   = 'Battery_ok';
      event(end  ).sample = i;
    end

  case 'brainvision_vmrk'
    fid=fopen(filename,'rt');
    if fid==-1,
      error('cannot open BrainVision marker file')
    end
    line = [];
    while ischar(line) || isempty(line)
      line = fgetl(fid);
      if ~isempty(line) && ~(isnumeric(line) && line==-1)
        if strncmpi(line, 'Mk', 2)
          % this line contains a marker
          tok = tokenize(line, '=', 0);    % do not squeeze repetitions of the seperator
          if length(tok)~=2
            warning('skipping unexpected formatted line in BrainVision marker file');
          else
            % the line looks like "MkXXX=YYY", which is ok
            % the interesting part now is in the YYY, i.e. the second token
            tok = tokenize(tok{2}, ',', 0);    % do not squeeze repetitions of the seperator
            if isempty(tok{1})
              tok{1}  = [];
            end
            if isempty(tok{2})
              tok{2}  = [];
            end
            event(end+1).type     = tok{1};
            event(end  ).value    = tok{2};
            event(end  ).sample   = str2num(tok{3});
            event(end  ).duration = str2num(tok{4});
          end
        end
      end
    end
    fclose(fid);

  case 'ced_son'
    % check that the required low-level toolbox is available
    hastoolbox('neuroshare', 1);
    orig = read_ced_son(filename,'readevents','yes');
    event = struct('type',     {orig.events.type},...
      'sample',   {orig.events.sample},...
      'value',    {orig.events.value},...
      'offset',   {orig.events.offset},...
      'duration', {orig.events.duration});

  case {'ctf_ds', 'ctf_meg4', 'ctf_res4'}
    % obtain the dataset name
    if filetype(filename, 'ctf_meg4') ||  filetype(filename, 'ctf_res4')
      filename = fileparts(filename);
    end
    [path, name, ext] = fileparts(filename);
    headerfile = fullfile(path, [name ext], [name '.res4']);
    datafile   = fullfile(path, [name ext], [name '.meg4']);
    classfile  = fullfile(path, [name ext], 'ClassFile.cls');
    markerfile = fullfile(path, [name ext], 'MarkerFile.mrk');

    % read the header, required to determine the stimulus channels and trial specification
    if isempty(hdr)
      hdr = read_header(headerfile);
    end

    % read the trigger codes from the STIM channel, usefull for (pseudo) continuous data
    % this splits the trigger channel into the lowers and highest 16 bits,
    % corresponding with the front and back panel of the electronics cabinet at the Donders Centre
    [backpanel, frontpanel] = read_ctf_trigger(filename);
    for i=find(backpanel(:)')
      event(end+1).type   = 'backpanel trigger';
      event(end  ).sample = i;
      event(end  ).value  = backpanel(i);
    end
    for i=find(frontpanel(:)')
      event(end+1).type   = 'frontpanel trigger';
      event(end  ).sample = i;
      event(end  ).value  = frontpanel(i);
    end

    % specify the range to search for triggers, default is the complete file
    if ~isempty(flt_minsample)
      begsample = flt_minsample;
    else
      begsample = 1;
    end
    if ~isempty(flt_maxsample)
      endsample = flt_maxsample;
    else
      endsample = hdr.nSamples*hdr.nTrials;
    end

    % determine the trigger channels from the header
    if isfield(hdr, 'orig') && isfield(hdr.orig, 'sensType')
      for i=find(hdr.orig.sensType(:)'==11)
        % read the trigger channel as raw data, can safely assume that it is continuous
        trig = read_data(datafile, hdr, begsample, endsample, i, 1);
        % correct for reading it as signed integer, whereas it should be an unsigned int
        trig(find(trig<0)) = trig(find(trig<0)) + 2^32;
        % convert the trigger into an event with a value at a specific sample
        for j=find(diff([0 trig(:)'])>0)
          event(end+1).type   = hdr.label{i};
          event(end  ).sample = j + begsample - 1;
          event(end  ).value  = trig(j);
        end
      end
    end

    % make an event for each trial as defined in the header
    for i=1:hdr.nTrials
      event(end+1).type     = 'trial';
      event(end  ).sample   = (i-1)*hdr.nSamples + 1;
      event(end  ).offset   = -hdr.nSamplesPre;
      event(end  ).duration =  hdr.nSamples;
      event(end  ).value    =  [];
    end

    % read the classification file and make an event for each classified trial
    [condNumbers,condLabels] = read_ctf_cls(classfile);
    if ~isempty(condNumbers)
      Ncond = length(condLabels);
      for i=1:Ncond
        for j=1:length(condNumbers{i})
          event(end+1).type     = 'classification';
          event(end  ).value    = condLabels{i};
          event(end  ).sample   = (condNumbers{i}(j)-1)*hdr.nSamples + 1;
          event(end  ).offset   = -hdr.nSamplesPre;
          event(end  ).duration =  hdr.nSamples;
        end
      end
    end

    if exist(markerfile)
      % read the marker file and make an event for each marker
      % this depends on the readmarkerfile function that I got from Tom Holroyd
      % I have not tested this myself extensively, since at the FCDC we
      % don't use the marker files
      mrk = readmarkerfile(filename);
      for i=1:mrk.number_markers
        for j=1:mrk.number_samples(i)
          % determine the location of the marker, expressed in samples
          trialnum = mrk.trial_times{i}(j,1);
          synctime = mrk.trial_times{i}(j,2);
          begsample = (trialnum-1)*hdr.nSamples + 1;    % of the trial, relative to the start of the datafile
          endsample = (trialnum  )*hdr.nSamples;        % of the trial, relative to the start of the datafile
          offset    = round(synctime*hdr.Fs);           % this is the offset (in samples) relative to time t=0 for this trial
          offset    = offset + hdr.nSamplesPre;         % and time t=0 corrsponds with the nSamplesPre'th sample
          % store this marker as an event
          event(end+1).type    = mrk.marker_names{i};
          event(end ).value    = [];
          event(end ).sample   = begsample + offset;
          event(end ).duration = 0;
          event(end ).offset   = offset;
        end
      end
    end

  case 'ctf_shm'
    % contact Robert Oostenveld if you are interested in real-time acquisition on the CTF system
    % read the events from shared memory
    event = read_shm_event(filename, varargin{:});

  case 'eep_avr'
    % check that the required low-level toolbox is available
    hastoolbox('eeprobe', 1);
    % the headerfile and datafile are the same
    if isempty(hdr)
      hdr = read_header(filename);
    end
    event(end+1).type     = 'average';
    event(end  ).sample   = 1;
    event(end  ).duration = hdr.nSamples;
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).value    = [];

  case 'eep_cnt'
    % check that the required low-level toolbox is available
    hastoolbox('eeprobe', 1);
    % try to read external trigger file in EEP format
    trgfile = [filename(1:(end-3)), 'trg'];
    if exist(trgfile, 'file')
      if isempty(hdr)
        hdr = read_header(filename);
      end
      tmp = read_eep_trg(trgfile);
      % translate the EEProbe trigger codes to events
      for i=1:length(tmp)
        event(i).type     = 'trigger';
        event(i).sample   = round((tmp(i).time/1000) * hdr.Fs) + 1;    % convert from ms to samples
        event(i).value    = tmp(i).code;
        event(i).offset   = 0;
        event(i).duration = 0;
      end
    else
      warning('no triggerfile was found');
    end

  case 'egi_egis'
    if isempty(hdr)
      hdr = read_header(filename);
    end
    fhdr   = hdr.orig.fhdr;
    chdr   = hdr.orig.chdr;
    ename  = hdr.orig.ename;
    cnames = hdr.orig.cnames;
    fcom   = hdr.orig.fcom;
    ftext  = hdr.orig.ftext;
    eventCount=0;
    for cell=1:fhdr(18)
      for trial=1:chdr(cell,2)
        eventCount=eventCount+1;
        event(eventCount).type     = 'trial';
        event(eventCount).sample   = (eventCount-1)*hdr.nSamples + 1;
        event(eventCount).offset   = -hdr.nSamplesPre;
        event(eventCount).duration =  hdr.nSamples;
        event(eventCount).value    =  cnames{cell};
      end
    end

  case 'egi_egia'
    if isempty(hdr)
      hdr = read_header(filename);
    end
    fhdr   = hdr.orig.fhdr;
    chdr   = hdr.orig.chdr;
    ename  = hdr.orig.ename;
    cnames = hdr.orig.cnames;
    fcom   = hdr.orig.fcom;
    ftext  = hdr.orig.ftext;
    eventCount=0;
    for cell=1:fhdr(18)
      for subject=1:chdr(cell,2)
        eventCount=eventCount+1;
        event(eventCount).type     = 'trial';
        event(eventCount).sample   = (eventCount-1)*hdr.nSamples + 1;
        event(eventCount).offset   = -hdr.nSamplesPre;
        event(eventCount).duration =  hdr.nSamples;
        event(eventCount).value    =  ['S' num2str(subject) cnames{cell}];
      end
    end

  case 'fcdc_matbin'
    % this is multiplexed data in a *.bin file, accompanied by a matlab file containing the header and event
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file '.mat']);
    % read the events from the Matlab file
    tmp   = load(filename, 'event');
    event = tmp.event;

  case 'fcdc_fifo'
    fifo = filetype_check_uri(filename);
    fid = fopen(fifo, 'r');
    msg = fread(fid, inf, 'char=>char');
    event = msg2struct(msg);
    % convert the message into event structure
    fclose(fid);

  case 'fcdc_tcpsocket'
    % TCP network socket
    [host, port] = filetype_check_uri(filename);
    con = pnet('tcpconnect', host, port);
    if con<0
      error('problem opening network connection');
    end
    % FIXME implementation should be finished
    keyboard

  case 'fcdc_serial'
    % serial port on windows or linux platform
    [port, opt] = filetype_check_uri(filename);
    % determine whether any serial port objects are already associated with the
    % target serial port
    s = [];
    temp = instrfind;
    if isa(temp,'instrument')
      % find all serial ports
      i1 = strcmp(lower({temp(:).Type}),'serial');
      if any(i1)
        % find all serial ports whose name matches that of the specified port
        i2 = strmatch(lower(port),lower({temp(find(i1)).Name}));
        % set s to the (first) matching port if present (and open if necessary)
        if ~isempty(i2)
          s = temp(i2(1));
          if ~strcmp(s.Status,'open'), fopen(s); end;
        end
      end
    end
    % create, configure a serial port object if necessary and open the port
    if ~isa(s,'serial')
      s = serial(port);
      if ~isempty(opt) && iscell(opt), s = set(s,opt); end;
      fopen(s);
    end
    % try to read a message from the serial port
    msg = [];
    % FIXME: this currently assumes that all messages are terminated by the
    % "newline" character (ascii character 10)
    try
      msg = fscanf(s,'%s\n');
    end;
    % convert message to event structure
    event = msg2struct(msg);

  case 'fcdc_mysql'
    % read from a MySQL server listening somewhere else on the network
    db_open(filename);
    if db_blob
      event = db_select_blob('fieldtrip.event', 'msg');
    else
      event = db_select('fieldtrip.event', {'type', 'value', 'sample', 'offset', 'duration'});
    end

  case 'matlab'
    % read the events from a normal Matlab file
    tmp   = load(filename, 'event');
    event = tmp.event;

  case {'mpi_ds', 'mpi_dap'}
    if isempty(hdr)
      hdr = read_header(filename);
    end
    % determine the DAP files that compromise this dataset
    if isdir(filename)
      ls = dir(filename);
      dapfile = {};
      for i=1:length(ls)
        if ~isempty(regexp(ls(i).name, '.dap$'))
          dapfile{end+1} = fullfile(filename, ls(i).name);
        end
      end
      dapfile = sort(dapfile);
    elseif iscell(filename)
      dapfile = filename;
    else
      dapfile = {filename};
    end
    % assume that each DAP file is accompanied by a dat file
    % read the trigger values from the separate dat files
    trg = [];
    for i=1:length(dapfile)
      datfile = [dapfile{i}(1:(end-4)) '.dat'];
      trg = cat(1, trg, textread(datfile, '', 'headerlines', 1));
    end
    % construct a event structure, one 'trialcode' event per trial
    for i=1:length(trg)
      event(i).type     = 'trialcode';            % string
      event(i).sample   = (i-1)*hdr.nSamples + 1; % expressed in samples, first sample of file is 1
      event(i).value    = trg(i);                 % number or string
      event(i).offset   = 0;                      % expressed in samples
      event(i).duration = hdr.nSamples;           % expressed in samples
    end

  case 'neuromag_fif'
    % check that the required low-level toolbox is available
    hastoolbox('meg-pd', 1);
    if isempty(hdr)
      hdr = read_header(filename);
    end
    % add the trials to the event structure
    for i=1:hdr.nTrials
      event(end+1).type     = 'trial';
      event(end  ).sample   = (i-1)*hdr.nSamples + 1;
      event(end  ).value    = [];
      event(end  ).offset   = -hdr.nSamplesPre;
      event(end  ).duration =  hdr.nSamples;
    end
    % add triggers based on the binary trigger channel, this is based the following email
    %
    % On 14 Sep 2004, at 16:33, Lauri Parkkonen wrote:
    %   Anke's file is probably quite recent one. It should have normal binary
    %   coding on STI014 and its "analog" counterparts on STI1 ... STI6,
    %   swinging between 0 ... +5 volts. These latter channels are virtual,
    %   i.e., they are generated from STI014 for backwards compatibility. STI1
    %   is the LSB of STI014 etc. For all the new stuff, I would recommend
    %   using STI014 as that's the way the triggers will be stored in the future
    %   (maybe the channel name is going to change to something more reasonable
    %   like TRIG). Unfortunately, some older files from the 306-channel system
    %   have the STI014 coded in a different way - the two 8-bit halves code
    %   the input and output triggers separately.
    triglab  = {'STI 014', 'STI 015', 'STI 016'};
    trigindx = match_str(hdr.label, triglab);
    if length(trigindx)>0
      % pad with the last sample of the previous block, start with zeros
      pad = zeros(length(trigindx),1);
      for i=1:hdr.nTrials
        begsample = (i-1)*hdr.nSamples + 1;
        endsample = (i  )*hdr.nSamples;
        raw = read_data(filename, hdr, begsample, endsample, trigindx);
        % detect the flank of the trigger
        flank = ((diff([raw pad], [], 2) .* raw)~=0);
        % for each flank in the TRIG channels, add an event
        for j=1:length(trigindx)
          for k=find(flank(j,:))
            event(end+1).type     = 'trigger';
            event(end  ).sample   = begsample + k;
            event(end  ).value    = raw(j,k);
            event(end  ).offset   = [];
            event(end  ).duration = [];
          end
        end
        % remember the last sample to ensure continuity with the next block
        pad = raw(:,end);
      end
    end % if length(trigindx)
    % look for the analog triggers only if no binary trigger channel is present
    if isempty(trigindx)
      % add the triggers to the event structure based on trigger channels with the name "STI xxx"
      % this is not a very efficient implementation, since it has to read all data
      % furthermore, there are some issues with noise on these analog trigger channels
      stimindx = [];
      stimlab  = {};
      for i=1:length(hdr.label)
        if all(hdr.label{i}(1:3) == 'STI')
          stimindx(end+1) = i;
          stimlab{end+1}  = hdr.label{i};
        end
      end
      if length(stimindx)>0
        % pad with the last sample of the previous block, start with zeros
        pad = zeros(length(stimindx),1);
        for i=1:hdr.nTrials
          begsample = (i-1)*hdr.nSamples + 1;
          endsample = (i  )*hdr.nSamples;
          raw = read_data(filename, hdr, begsample, endsample, stimindx);
          % detect the flank of the trigger
          flank = ((diff([raw pad], [], 2) .* raw)~=0) & (raw>=5);
          % Note: the original code read
          %   flank = (diff([raw pad], [], 2) .* raw)~=0;
          % but according to Joachim, real events always have triggers > 5
          % for each flank in the STIM channels, add an event
          for j=1:length(stimindx)
            for k=find(flank(j,:))
              event(end+1).type     = stimlab{j};
              event(end  ).sample   = begsample + k;
              event(end  ).value    = raw(j,k);
              event(end  ).offset   = [];
              event(end  ).duration = [];
            end
          end
          % remember the last sample to ensure continuity with the next block
          pad = raw(:,end);
        end
      end % if length(stimindx)
    end

  case {'neuralynx_ttl' 'neuralynx_bin' 'neuralynx_dma' 'neuralynx_sdma'}
    if isempty(hdr)
      hdr = read_header(filename);
    end
    
    % specify the range to search for triggers, default is the complete file
    if ~isempty(flt_minsample)
      begsample = flt_minsample;
    else
      begsample = 1;
    end
    if ~isempty(flt_maxsample)
      endsample = flt_maxsample;
    else
      endsample = hdr.nSamples*hdr.nTrials;
    end

    if strcmp(eventformat, 'neuralynx_dma')
      % read the Parallel_in channel from the DMA log file
      ttl = read_neuralynx_dma(filename, begsample, endsample, 'ttl');
    elseif strcmp(eventformat, 'neuralynx_sdma')
      % determine the seperate files with the trigger and timestamp information
      [p, f, x] = fileparts(filename);
      ttlfile = fullfile(filename, [f '.ttl.bin']);
      tslfile = fullfile(filename, [f '.tsl.bin']);
      tshfile = fullfile(filename, [f '.tsh.bin']);
      if ~exist(ttlfile) && ~exist(tslfile) && ~exist(tshfile)
        % perhaps it is an old splitted dma dataset?
        ttlfile = fullfile(filename, [f '.ttl']);
        tslfile = fullfile(filename, [f '.tsl']);
        tshfile = fullfile(filename, [f '.tsh']);
      end
      if ~exist(ttlfile) && ~exist(tslfile) && ~exist(tshfile)
        % these files must be present in a splitted dma dataset
        error('could not locate the individual ttl, tsl and tsh files');
      end
      % read the trigger values from the seperate file
      ttl = read_neuralynx_bin(ttlfile, begsample, endsample);
    elseif strcmp(eventformat, 'neuralynx_ttl')
      % determine the optional files with timestamp information
      tslfile = [filename(1:(end-4)) '.tsl'];
      tshfile = [filename(1:(end-4)) '.tsh'];
      % read the triggers from a seperate *.ttl file
      ttl = read_neuralynx_ttl(filename, begsample, endsample);
    elseif strcmp(eventformat, 'neuralynx_bin')
      % determine the optional files with timestamp information
      tslfile = [filename(1:(end-8)) '.tsl.bin'];
      tshfile = [filename(1:(end-8)) '.tsh.bin'];
      % read the triggers from a seperate *.ttl.bin file
      ttl = read_neuralynx_bin(filename, begsample, endsample);
    end

    ttl = int32(ttl / (2^16));    % parallel port provides int32, but word resolution is int16. Shift the bits and typecast to signed integer.
    d1  = (diff(ttl)~=0);         % determine the flanks, which can be multiple samples long (this looses one sample)
    d2  = (diff(d1)==1);          % determine the onset of the flanks (this looses one sample)
    smp = find(d2)+2;             % find the onset of the flanks, add the two samples again
    val = ttl(smp+5);             % look some samples further for the trigger value, to avoid the flank
    clear d1 d2 ttl
    ind = find(val~=0);           % look for triggers tith a non-zero value, this avoids downgoing flanks going to zero
    smp = smp(ind);               % discard triggers with a value of zero
    val = val(ind);               % discard triggers with a value of zero

    if ~isempty(smp)
      % try reading the timestamps
      if strcmp(eventformat, 'neuralynx_dma')
        tsl = read_neuralynx_dma(filename, 1, max(smp), 'tsl');
        tsl = tsl(smp);
        tsh = read_neuralynx_dma(filename, 1, max(smp), 'tsh');
        tsh = tsh(smp);
        ts  = timestamp_neuralynx(tsl, tsh);
      elseif exist(tslfile) && exist(tshfile)
        tsl = read_neuralynx_dma(filename, 1, max(smp), 'tsl');
        tsl = tsl(smp);
        tsh = read_neuralynx_dma(filename, 1, max(smp), 'tsh');
        tsh = tsh(smp);
        ts  = timestamp_neuralynx(tsl, tsh);
      else
        ts = [];
      end

      % reformat the values as cell array, since the struct function can work with those
      type      = repmat({'trigger'},size(smp));
      value     = num2cell(val);
      sample    = num2cell(smp + begsample - 1);
      duration  = repmat({[]},size(smp));
      offset    = repmat({[]},size(smp));
      if ~isempty(ts)
        value     = num2cell(ts);
      else
        timestamp  = repmat({[]},size(smp));
      end
      % convert it into a structure array, this can be done in one go
      event = struct('type', type, 'value', value, 'sample', sample, 'timestamp', timestamp, 'offset', offset, 'duration', duration);
      clear type value sample timestamp offset duration
    end

    if (strcmp(filetype, 'neuralynx_bin') || strcmp(filetype, 'neuralynx_ttl')) && isfield(hdr, 'FirstTimeStamp')
      % the header was obtained from an external dataset which could be at a different sampling rate
      % use the timestamps to redetermine the sample numbers
      fprintf('using sample number of the downsampled file to reposition the TTL events\n');
      % convert the timestamps into samples, keeping in mind the FirstTimeStamp and TimeStampPerSample
      smp = round(double(ts - uint64(hdr.FirstTimeStamp))./hdr.TimeStampPerSample + 1);
      for i=1:length(event)
        % update the sample number
        event(i).sample = smp(i);
      end
    end

  case 'neuralynx_ds'
    % read the header of the dataset
    if isempty(hdr)
      hdr = read_header(filename);
    end
    % the event file is contained in the dataset directory
    if     exist(fullfile(filename, 'Events.Nev'))
      filename = fullfile(filename, 'Events.Nev');
    elseif exist(fullfile(filename, 'Events.nev'))
      filename = fullfile(filename, 'Events.nev');
    elseif exist(fullfile(filename, 'events.Nev'))
      filename = fullfile(filename, 'events.Nev');
    elseif exist(fullfile(filename, 'events.nev'))
      filename = fullfile(filename, 'events.nev');
    end
    % read the events, apply filter is applicable
    nev = read_neuralynx_nev(filename, 'type', flt_type, 'value', flt_value, 'mintimestamp', flt_mintimestamp, 'maxtimestamp', flt_maxtimestamp, 'minnumber', flt_minnumber, 'maxnumber', flt_maxnumber);

    % now get the values as cell array, since the struct function can work with those
    value     = {nev.TTLValue};
    timestamp = {nev.TimeStamp};
    number    = {nev.EventNumber};
    type      = repmat({'trigger'},size(value));
    duration  = repmat({[]},size(value));
    offset    = repmat({[]},size(value));
    sample    = num2cell(round(double(cell2mat(timestamp) - hdr.FirstTimeStamp)/hdr.TimeStampPerSample + 1));
    % convert it into a structure array
    event = struct('type', type, 'value', value, 'sample', sample, 'timestamp', timestamp, 'duration', duration, 'offset', offset, 'number', number);

  case 'neuralynx_cds'
    % this is a combined Neuralynx dataset with seperate subdirectories for the LFP, MUA and spike channels
    dirlist   = dir(filename);
    %haslfp   = any(filetype_check_extension({dirlist.name}, 'lfp'));
    %hasmua   = any(filetype_check_extension({dirlist.name}, 'mua'));
    %hasspike = any(filetype_check_extension({dirlist.name}, 'spike'));
    %hastsl   = any(filetype_check_extension({dirlist.name}, 'tsl'));   % seperate file with original TimeStampLow
    %hastsh   = any(filetype_check_extension({dirlist.name}, 'tsh'));   % seperate file with original TimeStampHi
    hasttl    = any(filetype_check_extension({dirlist.name}, 'ttl'));   % seperate file with original Parallel_in
    hasnev    = any(filetype_check_extension({dirlist.name}, 'nev'));   % original Events.Nev file
    hasmat    = 0;
    if hasttl
      eventfile = fullfile(filename, dirlist(find(filetype_check_extension({dirlist.name}, 'ttl'))).name);
      % read the header from the combined dataset
      if isempty(hdr)
        hdr = read_header(filename);
      end
      % read the events from the *.ttl file
      event = read_event(eventfile);
      % convert the sample numbers from the dma or ttl file to the downsampled dataset
      % assume that the *.ttl file is sampled at 32556Hz and is aligned with the rest of the data
      for i=1:length(event)
        event(i).sample = round((event(i).sample-1) * hdr.Fs/32556 + 1);
      end
      % elseif hasnev
      % FIXME, do something here
      % elseif hasmat
      % FIXME, do something here
    else
      error('no event file found');
    end

%   The sample number is missingin the code below, since it is not available
%   without looking in the continuously sampled data files. Therefore
%   sorting the events (later in this function) based on the sample number
%   fails and no events can be returned.
%
%   case 'neuralynx_nev'
%     [nev] = read_neuralynx_nev(filename);
%     % select only the events with a TTL value
%     ttl = [nev.TTLValue];
%     sel = find(ttl~=0);
%     % now get the values as cell array, since teh struct function can work with those
%     value     = {nev(sel).TTLValue};
%     timestamp = {nev(sel).TimeStamp};
%     event = struct('value', value, 'timestamp', timestamp);
%     for i=1:length(event)
%       % assign the other fixed elements
%       event(i).type     = 'trigger';
%       event(i).offset   = [];
%       event(i).duration = [];
%       event(i).sample   = [];
%     end

  case 'nexstim_nxe'
    event = read_nexstim_event(filename);

  case 'ns_avg'
    if isempty(hdr)
      hdr = read_header(filename);
    end
    event(end+1).type     = 'average';
    event(end  ).sample   = 1;
    event(end  ).duration = hdr.nSamples;
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).value    = [];

  case 'ns_cnt'
    % read the header, the original header includes the event table
    if isempty(hdr)
      hdr = read_header(filename);
    end
    % translate the event table into known FieldTrip event types
    for i=1:hdr.orig.nevent
      event(i).type     = 'trigger';
      event(i).sample   = hdr.orig.event.frame(i);
      event(i).value    = hdr.orig.event.stimtype(i);
      event(i).offset   = 0;
      event(i).duration = 0;
    end

  case 'ns_eeg'
    if isempty(hdr)
      hdr = read_header(filename);
    end
    for i=1:hdr.nTrials
      event(end+1).type     = 'trial';
      event(end  ).sample   = (i-1)*hdr.nSamples + 1;
      event(end  ).value    = [];
      event(end  ).offset   = -hdr.nSamplesPre;
      event(end  ).duration =  hdr.nSamples;
      % read the data to determine manually accepted/rejected trials
      tmp = read_ns_eeg(filename, i);
      if tmp.sweep.accept
        event(end).value = 'accept';
      else
        event(end).value = 'reject';
      end
    end

  case 'plexon_nex'
    event = read_nex_event(filename);

  case 'yokogawa_ave'
    % check that the required low-level toolbox is available
    hastoolbox('yokogawa', 1);
    if isempty(hdr)
      hdr = read_header(filename);
    end
    event(end+1).type     = 'average';
    event(end  ).sample   = 1;
    event(end  ).duration = hdr.nSamples;
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).value    = [];

  case 'yokogawa_con'
    % check that the required low-level toolbox is available
    % hastoolbox('yokogawa', 1);
    error('events still need to be implemented for the yokogawa_con format');

  case 'yokogawa_raw'
    % check that the required low-level toolbox is available
    hastoolbox('yokogawa', 1);
    % read the trigger id from all trials
    value = GetMeg160TriggerEventM(filename);
    % create a "trial" event for each trial and assign it the corresponding trigger value
    for i=1:hdr.nTrials
      event(end+1).type     = 'trial';
      event(end  ).sample   = (i-1)*hdr.nSamples + 1;
      event(end  ).offset   = -hdr.nSamplesPre;
      event(end  ).duration =  hdr.nSamples;
      event(end  ).value    = value(i);
    end

  otherwise
    error('unsupported event format');
end

if ~isempty(event)
  % make sure that all required elements are present
  if ~isfield(event, 'type'),     error('type field not defined for each event');     end
  if ~isfield(event, 'sample'),   error('sample field not defined for each event');   end
  if ~isfield(event, 'value'),    for i=1:length(event), event(i).value = [];    end; end
  if ~isfield(event, 'offset'),   for i=1:length(event), event(i).offset = [];   end; end
  if ~isfield(event, 'duration'), for i=1:length(event), event(i).duration = []; end; end
end

if ~isempty(event)
  % sort the events on the sample on which they occur
  % this has the side effect that events without a sample number are discarded
  [dum, indx] = sort([event.sample]);
  event = event(indx);
else
  warning(sprintf('no events found in %s', filename));
end

% apply the optional filters
event = filter_event(event, varargin{:});

function c = catevent(a, b);
if isempty(a)
  c = b(:);
elseif isempty(b)
  c = a(:);
else
  c = a(:);
  for i=1:numel(b)
    c(end+1).type     = b(i).type;
    c(end  ).sample   = b(i).sample;
    if isfield(b, 'timestamp')
      c(end  ).timestamp= b(i).timestamp;
    end
    c(end  ).value    = b(i).value;
    c(end  ).offset   = b(i).offset;
    c(end  ).duration = b(i).duration;
  end
end


