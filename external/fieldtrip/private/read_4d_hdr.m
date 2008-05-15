function [header] = read_4d_hdr(datafile, configfile)

%hdr=READ_4D_HDR(filename)
%Collects the required Fieldtrip header data from the data file 'filename'
%and the associated 'config' file for that data.

%Gavin Paterson 28/01/08 
%g.paterson@psy.gla.ac.uk
%
%Adapted from the MSI>>Matlab code written by:
%eugene.kronberg@uchsc.edu


%read header
if nargin ~= 2
    error('Wrong number of input arguments');
end


%always big endian  
fid = fopen(datafile, 'r', 'b');

if fid == -1
    error('Cannot open file %s', datafile);
end

fseek(fid, 0, 'eof');
header_end = ftell(fid);
%last 8 bytes of the pdf is header offset
fseek(fid, -8, 'eof');
header_offset = fread(fid,1,'uint64');

%first byte of the header
fseek(fid, header_offset, 'bof');

% read header data
align_file_pointer(fid)
header.header_data.FileType  = fread(fid, 1, 'uint16=>uint16');
file_type                    = char(fread(fid, 5, 'uchar'))';
header.header_data.file_type = file_type(file_type>0);
fseek(fid, 1, 'cof');
format = fread(fid, 1, 'int16=>int16');
switch format
    case 1
        header.header_data.Format = 'SHORT';
    case 2
        header.header_data.Format = 'LONG';
    case 3
        header.header_data.Format = 'FLOAT';
    case 4
        header.header_data.Format ='DOUBLE';
end
header.header_data.acq_mode           = fread(fid, 1,  'uint16=>uint16');
header.header_data.TotalEpochs        = fread(fid, 1,  'uint32=>double');
header.header_data.input_epochs       = fread(fid, 1,  'uint32=>uint32');
header.header_data.TotalEvents        = fread(fid, 1,  'uint32=>uint32');
header.header_data.total_fixed_events = fread(fid, 1,  'uint32=>uint32');
header.header_data.SamplePeriod       = fread(fid, 1,  'float32=>float64');
header.header_data.SampleFrequency    = 1/header.header_data.SamplePeriod;
xaxis_label                           = char(fread(fid, 16, 'uchar'))';
header.header_data.xaxis_label        = xaxis_label(xaxis_label>0);
header.header_data.total_processes    = fread(fid, 1,  'uint32=>uint32');
header.header_data.TotalChannels      = fread(fid, 1,  'uint16=>double');
fseek(fid, 2, 'cof');
header.header_data.checksum           = fread(fid, 1,  'int32=>int32');
header.header_data.total_ed_classes   = fread(fid, 1,  'uint32=>uint32');
header.header_data.total_associated_files = fread(fid, 1, 'uint16=>uint16');
header.header_data.last_file_index    = fread(fid, 1,  'uint16=>uint16');
header.header_data.timestamp          = fread(fid, 1,  'uint32=>uint32');
header.header_data.reserved           = fread(fid, 20, 'uchar')';
fseek(fid, 4, 'cof');
 
%read epoch_data
for epoch = 1:header.header_data.TotalEpochs;
    align_file_pointer(fid)

    header.epoch_data(epoch).pts_in_epoch     = fread(fid, 1,  'uint32=>uint32');
    header.epoch_data(epoch).epoch_duration   = fread(fid, 1,  'float32=>float32');
    header.epoch_data(epoch).expected_iti     = fread(fid, 1,  'float32=>float32');
    header.epoch_data(epoch).actual_iti       = fread(fid, 1,  'float32=>float32');
    header.epoch_data(epoch).total_var_events = fread(fid, 1,  'uint32=>uint32');
    header.epoch_data(epoch).checksum         = fread(fid, 1,  'int32=>int32');
    header.epoch_data(epoch).epoch_timestamp  = fread(fid, 1,  'int32=>int32');
    header.epoch_data(epoch).reserved         = fread(fid, 28, 'uchar')';
    header.header_data.SlicesPerEpoch         = double(header.epoch_data(1).pts_in_epoch);

    %read event data (var_events)
    for event = 1:header.epoch_data(epoch).total_var_events
        align_file_pointer(fid)

        event_name = char(fread(fid, 16, 'uchar'))';
        header.epoch_data(epoch).var_event{event}.event_name  = event_name(event_name>0);
        header.epoch_data(epoch).var_event{event}.start_lat   = fread(fid, 1,  'float32=>float32');
        header.epoch_data(epoch).var_event{event}.end_lat     = fread(fid, 1,  'float32=>float32');
        header.epoch_data(epoch).var_event{event}.step_size   = fread(fid, 1,  'float32=>float32');
        header.epoch_data(epoch).var_event{event}.fixed_event = fread(fid, 1,  'uint16=>uint16');
        fseek(fid, 2, 'cof');
        header.epoch_data(epoch).var_event{event}.checksum    = fread(fid, 1,  'int32=>int32');
        header.epoch_data(epoch).var_event{event}.reserved    = fread(fid, 32, 'uchar')';
        fseek(fid, 4, 'cof');
    end
end

%read  channel ref data
for channel = 1:header.header_data.TotalChannels
    align_file_pointer(fid)

    chan_label                                 = (fread(fid, 16, 'uchar=>char'))';
    header.channel_data(channel).chan_label    = chan_label(chan_label>0);
    header.channel_data(channel).chan_no       = fread(fid, 1, 'uint16=>uint16');
    header.channel_data(channel).attributes    = fread(fid, 1, 'uint16=>uint16');
    header.channel_data(channel).scale         = fread(fid, 1, 'float32=>float32');
    yaxis_label                                = char(fread(fid, 16, 'uchar=>char'))';
    header.channel_data(channel).yaxis_label   = yaxis_label(yaxis_label>0);
    header.channel_data(channel).valid_min_max = fread(fid, 1, 'uint16=>uint16');
    fseek(fid, 6, 'cof');
    header.channel_data(channel).ymin          = fread(fid, 1,  'float64');
    header.channel_data(channel).ymax          = fread(fid, 1,  'float64');
    header.channel_data(channel).index         = fread(fid, 1,  'uint32=>uint32');
    header.channel_data(channel).checksum      = fread(fid, 1,  'int32=>int32');
    header.channel_data(channel).whatisit      = char(fread(fid, 4, 'uchar=>char'))';
    header.channel_data(channel).reserved      = fread(fid, 28, 'uchar')';
end

%read event data
for event = 1:header.header_data.total_fixed_events
    align_file_pointer(fid)
    event_name                           = char(fread(fid, 16, 'uchar'))';
    header.event_data(event).event_name  = event_name(event_name>0);
    header.event_data(event).start_lat   = fread(fid, 1,  'float32=>float32');
    header.event_data(event).end_lat     = fread(fid, 1,  'float32=>float32');
    header.event_data(event).step_size   = fread(fid, 1,  'float32=>float32');
    header.event_data(event).fixed_event = fread(fid, 1,  'uint16=>uint16');
    fseek(fid, 2, 'cof');
    header.event_data(event).checksum    = fread(fid, 1,  'int32=>int32');
    header.event_data(event).reserved    = fread(fid, 32, 'uchar')';
    fseek(fid, 4, 'cof');
end
header.header_data.FirstLatency = double(header.event_data(1).start_lat);

%experimental: read process information
for np = 1:header.header_data.total_processes
  align_file_pointer(fid)
  nbytes                          = fread(fid, 1,  'uint32=>uint32');
  fp                              = ftell(fid);
  header.process(np).hdr.nbytes   = nbytes;
  type                            = char(fread(fid, 20, 'uchar'))';
  header.process(np).hdr.type     = type(type>0);
  header.process(np).hdr.checksum = fread(fid, 1,  'int32=>int32'); 
  user                            = char(fread(fid, 32, 'uchar'))';
  header.process(np).user         = user(user>0);
  header.process(np).timestamp    = fread(fid, 1,  'uint32=>uint32');
  fname                           = char(fread(fid, 32, 'uchar'))';
  header.process(np).filename     = fname(fname>0);
  fseek(fid, 28*8, 'cof'); %dont know
  header.process(np).totalsteps   = fread(fid, 1,  'uint32=>uint32');
  header.process(np).checksum     = fread(fid, 1,  'int32=>int32');
  header.process(np).reserved     = fread(fid, 32, 'uchar')'; 
  for ns = 1:header.process(np).totalsteps
    header.process(np).step(ns).hdr.nbytes    = fread(fid, 1, 'uint32=>uint32');
    type                                      = char(fread(fid, 20, 'uchar'))';
    header.process(np).step(ns).hdr.type      = type(type>0); %dont know how to interpret the first two
    header.process(np).step(ns).hdr.checksum  = fread(fid, 1, 'int32=>int32');
    header.process(np).step(ns).userblocksize = fread(fid, 1, 'int32=>int32');
    fseek(fid, 32, 'cof'); %needed until next step FIXME make more robust, the total number of read bytes
    %should be equal to the nbytes computed earlier on
    if strcmp(header.process(np).step(ns).hdr.type, 'PDF_Weight_Table'),
      warning('reading in weight table: no warranty that this is correct. it seems to work for the Glasgow 248-magnetometer system. if you have some code yourself, and/or would like to test it on your own data, please contact Jan-Mathijs');
      tmpfp = ftell(fid);
      tmp   = fread(fid, 1, 'uint8');
      Nchan = fread(fid, 1, 'uint32');
      Nref  = fread(fid, 1, 'uint32');
      for k = 1:Nref
        name                                     = fread(fid, 17, 'uchar'); %strange, but true
        header.process(np).step(ns).RefChan{k,1} = char(name(name>0))';
      end
      fseek(fid, 152, 'cof');
      for k = 1:Nchan
        name                                     = fread(fid, 17, 'uchar');
	header.process(np).step(ns).Chan{k,1}   = char(name(name>0))';
      end
      %fseek(fid, 20, 'cof');
      %fseek(fid, 4216, 'cof');
      header.process(np).step(ns).stuff1  = fread(fid, 4236, 'uint8');
      name                                = fread(fid, 16, 'uchar');
      header.process(np).step(ns).Creator = char(name(name>0))';
      %some stuff I don't understand yet
      %fseek(fid, 136, 'cof');
      header.process(np).step(ns).stuff2  = fread(fid, 136, 'uint8');
      %now something strange is going to happen: the weights are probably little-endian encoded.
      %here we go: check whether this applies to the whole PDF weight table
      fp = ftell(fid);
      fclose(fid);
      fid = fopen(datafile, 'r', 'l');
      fseek(fid, fp, 'bof');
      for k = 1:Nchan
        header.process(np).step(ns).Weights(k,:) = fread(fid, 23, 'float32=>float32')';
	fseek(fid, 36, 'cof');
      end
    else
    end    
  end
end
fclose(fid);

%end read header

%read config file
fid = fopen(configfile, 'r', 'b');

if fid == -1
    error('Cannot open config file');
end

header.config_data.version           = fread(fid, 1, 'uint16=>uint16');
site_name                            = char(fread(fid, 32, 'uchar'))';
header.config_data.site_name         = site_name(site_name>0);
dap_hostname                         = char(fread(fid, 16, 'uchar'))';
header.config_data.dap_hostname      = dap_hostname(dap_hostname>0);
header.config_data.sys_type          = fread(fid, 1, 'uint16=>uint16');
header.config_data.sys_options       = fread(fid, 1, 'uint32=>uint32');
header.config_data.supply_freq       = fread(fid, 1, 'uint16=>uint16');
header.config_data.total_chans       = fread(fid, 1, 'uint16=>uint16');
header.config_data.system_fixed_gain = fread(fid, 1, 'float32=>float32');
header.config_data.volts_per_bit     = fread(fid, 1, 'float32=>float32');
header.config_data.total_sensors     = fread(fid, 1, 'uint16=>uint16');
header.config_data.total_user_blocks = fread(fid, 1, 'uint16=>uint16');
header.config_data.next_derived_channel_number = fread(fid, 1, 'uint16=>uint16');
fseek(fid, 2, 'cof');
header.config_data.checksum          = fread(fid, 1, 'int32=>int32');
header.config_data.reserved          = fread(fid, 32, 'uchar=>uchar')';


header.config.Xfm = fread(fid, [4 4], 'double');

%user blocks
for ub = 1:header.config_data.total_user_blocks
    align_file_pointer(fid)
    header.user_block_data(ub).hdr.nbytes   = fread(fid, 1, 'uint32=>uint32');
    type                                    = char(fread(fid, 20, 'uchar'))';
    header.user_block_data(ub).hdr.type     = type(type>0);
    header.user_block_data(ub).hdr.checksum = fread(fid, 1, 'int32=>int32');
    user                                    = char(fread(fid, 32, 'uchar'))';
    header.user_block_data(ub).user         = user(user>0);
    header.user_block_data(ub).timestamp    = fread(fid, 1, 'uint32=>uint32');
    header.user_block_data(ub).user_space_size = fread(fid, 1, 'uint32=>uint32');
    header.user_block_data(ub).reserved     = fread(fid, 32, 'uchar=>uchar')';
    fseek(fid, 4, 'cof');
    user_space_size                         = double(header.user_block_data(ub).user_space_size);
    if strcmp(type(type>0), 'B_weights_used'), 
      %warning('reading in weight table: no warranty that this is correct');
      warning('reading in weight table: no warranty that this is correct. it seems to work for the Glasgow 248-magnetometer system. if you have some code yourself, and/or would like to test it on your own data, please contact Jan-Mathijs');
      tmpfp = ftell(fid);
      %read user_block_data weights
      %fseek(fid, 12, 'cof');
      %there is information in the 4th and 8th byte, these might be related to the settings?
      unknown1 = fread(fid, 1, 'uint32');
      unknown2 = fread(fid, 1, 'uint32');
      Nchan    = fread(fid, 1, 'uint32');
      Position = fread(fid, 32, 'uchar');
      header.user_block_data(ub).position = char(Position(Position>0))';
      fseek(fid, tmpfp+124, 'bof');
      Nanalog  = fread(fid, 1, 'uint32');
      Ndigital = fread(fid, 1, 'uint32');
      fseek(fid, tmpfp+204, 'bof');
      for k = 1:248
        Name     = fread(fid, 16, 'uchar');
        header.user_block_data(ub).channames{k,1} = char(Name(Name>0))';
      end
      for k = 1:Nanalog
        Name     = fread(fid, 16, 'uchar');
	header.user_block_data(ub).arefnames{k,1} = char(Name(Name>0))';
      end
      for k = 1:Ndigital
        Name     = fread(fid, 16, 'uchar');
	header.user_block_data(ub).drefnames{k,1} = char(Name(Name>0))';
      end

      header.user_block_data(ub).dweights = fread(fid, [Ndigital Nchan], 'single=>double')';
      header.user_block_data(ub).aweights = fread(fid, [Nanalog  Nchan],  'int16')'; 
      fseek(fid, tmpfp, 'bof');
    elseif strcmp(type(type>0), 'B_E_table_used'),
      %warning('reading in weight table: no warranty that this is correct');
      %tmpfp = ftell(fid);
      %fseek(fid, 4, 'cof'); %there's info here dont know how to interpret
      %Nx    = fread(fid, 1, 'uint32');
      %Nchan = fread(fid, 1, 'uint32');
      %type  = fread(fid, 32, 'uchar'); %don't know whether correct
      %header.user_block_data(ub).type = char(type(type>0))';
      %fseek(fid, 16, 'cof');
      %for k = 1:Nchan
      %  name                                 = fread(fid, 16, 'uchar');
      %  header.user_block_data(ub).name{k,1} = char(name(name>0))';
      %end
    end
    fseek(fid, user_space_size, 'cof');
end

%channels
for ch = 1:header.config_data.total_chans
    align_file_pointer(fid)
    name                                       = char(fread(fid, 16, 'uchar'))';
    header.config.channel_data(ch).name        = name(name>0);
    header.config.channel_data(ch).chan_no     = fread(fid, 1, 'uint16=>uint16');
    header.config.channel_data(ch).type        = fread(fid, 1, 'uint16=>uint16');
    header.config.channel_data(ch).sensor_no   = fread(fid, 1, 'int16=>int16');
    fseek(fid, 2, 'cof');
    header.config.channel_data(ch).gain        = fread(fid, 1, 'float32=>float32');
    header.config.channel_data(ch).units_per_bit = fread(fid, 1, 'float32=>float32');
    yaxis_label                                = char(fread(fid, 16, 'uchar'))';
    header.config.channel_data(ch).yaxis_label = yaxis_label(yaxis_label>0);
    header.config.channel_data(ch).aar_val     = fread(fid, 1, 'double');
    header.config.channel_data(ch).checksum    = fread(fid, 1, 'int32=>int32');
    header.config.channel_data(ch).reserved    = fread(fid, 32, 'uchar=>uchar')';
    fseek(fid, 4, 'cof');

    align_file_pointer(fid)
    header.config.channel_data(ch).device_data.hdr.size     = fread(fid, 1, 'uint32=>uint32');
    header.config.channel_data(ch).device_data.hdr.checksum = fread(fid, 1, 'int32=>int32');
    header.config.channel_data(ch).device_data.hdr.reserved = fread(fid, 32, 'uchar=>uchar')';

    switch header.config.channel_data(ch).type
        case {1, 3}%meg/ref

            header.config.channel_data(ch).device_data.inductance  = fread(fid, 1, 'float32=>float32');
            fseek(fid, 4, 'cof');
            header.config.channel_data(ch).device_data.Xfm         = fread(fid, [4 4], 'double');
            header.config.channel_data(ch).device_data.xform_flag  = fread(fid, 1, 'uint16=>uint16');
            header.config.channel_data(ch).device_data.total_loops = fread(fid, 1, 'uint16=>uint16');
            header.config.channel_data(ch).device_data.reserved    = fread(fid, 32, 'uchar=>uchar')';
            fseek(fid, 4, 'cof');

            for loop = 1:header.config.channel_data(ch).device_data.total_loops
                align_file_pointer(fid)
                header.config.channel_data(ch).device_data.loop_data(loop).position    = fread(fid, 3, 'double');
                header.config.channel_data(ch).device_data.loop_data(loop).direction   = fread(fid, 3, 'double');
                header.config.channel_data(ch).device_data.loop_data(loop).radius      = fread(fid, 1, 'double');
                header.config.channel_data(ch).device_data.loop_data(loop).wire_radius = fread(fid, 1, 'double');
                header.config.channel_data(ch).device_data.loop_data(loop).turns       = fread(fid, 1, 'uint16=>uint16');
                fseek(fid, 2, 'cof');
                header.config.channel_data(ch).device_data.loop_data(loop).checksum    = fread(fid, 1, 'int32=>int32');
                header.config.channel_data(ch).device_data.loop_data(loop).reserved    = fread(fid, 32, 'uchar=>uchar')';
            end
        case 2%eeg
            header.config.channel_data(ch).device_data.impedance       = fread(fid, 1, 'float32=>float32');
            fseek(fid, 4, 'cof');
            header.config.channel_data(ch).device_data.Xfm             = fread(fid, [4 4], 'double');
            header.config.channel_data(ch).device_data.reserved        = fread(fid, 32, 'uchar=>uchar')';
        case 4%external
            header.config.channel_data(ch).device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
            header.config.channel_data(ch).device_data.reserved        = fread(fid, 32, 'uchar=>uchar')';
            fseek(fid, 4, 'cof');
        case 5%TRIGGER
            header.config.channel_data(ch).device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
            header.config.channel_data(ch).device_data.reserved        = fread(fid, 32, 'uchar=>uchar')';
            fseek(fid, 4, 'cof');
        case 6%utility
            header.config.channel_data(ch).device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
            header.config.channel_data(ch).device_data.reserved        = fread(fid, 32, 'uchar=>uchar')';
            fseek(fid, 4, 'cof');
        case 7%derived
            header.config.channel_data(ch).device_data.user_space_size = fread(fid, 1, 'uint32=>uint32');
            header.config.channel_data(ch).device_data.reserved        = fread(fid, 32, 'uchar=>uchar')';
            fseek(fid, 4, 'cof');
        case 8%shorted
            header.config.channel_data(ch).device_data.reserved        = fread(fid, 32, 'uchar=>uchar')';
        otherwise
            error('Unknown device type: %d\n', header.config.channel_data(ch).type);
    end
end

fclose(fid);
%end read config file

header.header_data.FileDescriptor = 0; %no obvious field to take this from
header.header_data.Events         = 1;%no obvious field to take this from
header.header_data.EventCodes     = 0;%no obvious field to take this from

header.ChannelGain        = double([header.config.channel_data([header.channel_data.chan_no]).gain]');
header.ChannelUnitsPerBit = double([header.config.channel_data([header.channel_data.chan_no]).units_per_bit]');
header.Channel            = {header.config.channel_data([header.channel_data.chan_no]).name}';
header.Format             = header.header_data.Format;

function align_file_pointer(fid)
current_position = ftell(fid);
if mod(current_position, 8) ~= 0
    offset = 8 - mod(current_position,8);
    fseek(fid, offset, 'cof');
end
