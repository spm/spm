function write_data(filename, dat, varargin)

% WRITE_DATA exports electrophysiological data to a file
%
% Use as
%   write_data(filename, dat, ...)
%
% The specified filename can already contain the filename extention,
% but that is not required since it will be added automatically.
%
% Additional options should be specified in key-value pairs and can be
%   'header'         header structure, see READ_FCDC_HEADER
%   'dataformat'     string, see below
%   'chanindx'       1xN array
%   'subformat'      string, for the fcdc_ftc format
%
% The supported dataformats are
%   brainvision_eeg
%   matlab
%   riff_wave
%   fcdc_matbin
%   fcdc_mysql
%   plexon_nex
%   neuralynx_ncs
%   ctf_meg4       (partial and incomplete)
%
% See also READ_HEADER, READ_DATA, READ_EVEN, WRITE_EVENT

% Copyright (C) 2007-2008, Robert Oostenveld
%
% $Log: write_data.m,v $
% Revision 1.9  2008/01/31 20:05:05  roboos
% removed fcdc_ftc
% updated documentation
%
% Revision 1.8  2007/12/12 14:40:33  roboos
% added riff_wave using standard matlab function
%
% Revision 1.7  2007/12/06 17:09:06  roboos
% added fcdc_ftc
%
% Revision 1.6  2007/11/07 10:49:07  roboos
% cleaned up the reading and writing from/to mysql database, using db_xxx helper functions (see mysql directory)
%
% Revision 1.5  2007/11/05 17:02:45  roboos
% made first implementation for writing header and data to fcdc_mysql
%
% Revision 1.4  2007/10/01 13:49:11  roboos
% added skeleton implementation for ctf_meg4, not yet tested
%
% Revision 1.3  2007/06/13 08:14:58  roboos
% updated documentation
%
% Revision 1.2  2007/06/13 08:07:14  roboos
% added autodetection of dataformat
% added support for write_brainvision_eeg
%
% Revision 1.1  2007/06/13 06:44:13  roboos
% moved the content of the write_fcdc_data to the low-level write_data function
% updated the help
%

% for fcdc_mysql
global db_blob
if isempty(db_blob)
  db_blob = 0;
end

% get the options
dataformat    = keyval('dataformat',    varargin); if isempty(dataformat), dataformat = filetype(filename); end
subformat     = keyval('subformat',     varargin); % for fcdc_ftc
nbits         = keyval('nbits',         varargin); % for riff_wave
chanindx      = keyval('chanindx',      varargin);
hdr           = keyval('header',        varargin);

% determine the data size
[nchans, nsamples] = size(dat);

switch dataformat
  case 'ctf_meg4'  % wat als ctf_ds??
    % this is a skeleton implementation only and a lot of details still
    % need to be filled in. The implementation has not been tested yet.
    warning('this implementation has not yet been tested');
    
    if length(size(dat))<3
      dat = reshape(dat, [1 size(dat)]);
    end

    ntrldat   = size(dat,1);
    nsmpdat   = size(dat,3);
    nchandat  = size(dat,2);
    ntrlorig  = hdr.nTrials;
    nchanorig = hdr.nChans;
    nsmporig = hdr.nSamples;

    fid = fopen(filename, 'wb', 'ieee-be');
    buf = [77 69 71 52 49 67 80 0]; % 'MEG41CP', null-terminated
    fwrite(fid, buf, 'char');

    if ntrldat>ntrlorig
      error('');
    elseif nchandat>nchanorig
      error('');
    elseif ncsmpdat>nsmporig
      error('');
    end

    if isempty(chanindx) && length(hdr.label)==nchandat
      % vgl label met orig.label
      chanindx = 1:nchandat;
      % wat als nchandat~=nchanorig?
    end

    % wat als numbytes>1GB ?

    for i=1:ntrials
      datorig = zeros(nchanorig,nsamples)
      if i<=ntrldat
        datorig(chanindx,:) = dat(i,:,:); % padden met 0 als nsmpdat~=nsmporig? of continue maken en dan uitknippen?
        % wat als data>intmax?
        datorig = transpose(sparse(diag(1./hdr.orig.gainV)) * datorig);
        if any(datorig(:)>intmax('int32')) || any(datorig(:)<intmin('int32'))
          warning('data values were clipped to fit into 32 bit integer values');
          datorig(datorig>intmax('int32')) = intmax('int32');
          datorig(datorig<intmin('int32')) = intmin('int32');
        end
        datorig = int32(datorig);
      else
        datorig = int32(transpose(datorig));
      end
      fwrite(fid, datorig, 'int32');
    end

  case 'brainvision_eeg'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % combination of *.eeg and *.vhdr file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nchans~=hdr.nChans && length(chanindx)==nchans
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      hdr.label  = hdr.label(chanindx);
      hdr.nChans = length(chanindx);
    end
    % the header should at least contain the following fields
    %   hdr.label
    %   hdr.nChans
    %   hdr.Fs
    write_brainvision_eeg(filename, hdr, dat);

  case 'matlab'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plain matlab file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file '.mat']);
    save(filename, 'dat', 'hdr');

  case 'riff_wave'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     This writes data Y to a Windows WAVE file specified by the file name
    %     WAVEFILE, with a sample rate of FS Hz and with NBITS number of bits.
    %     NBITS must be 8, 16, 24, or 32.  For NBITS < 32, amplitude values
    %     outside the range [-1,+1] are clipped
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nchans~=hdr.nChans && length(chanindx)==nchans
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      hdr.label  = hdr.label(chanindx);
      hdr.nChans = length(chanindx);
    end
    if nchans~=1
      error('this format only supports single channel continuous data');
    end
    wavwrite(dat, hdr.Fs, nbits, filename);
    
  case 'fcdc_matbin'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multiplexed data in a *.bin file (ieee-le, 64 bit floating point values),
    % accompanied by a matlab V6 file containing the header
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.mat']);
    datafile   = fullfile(path, [file '.bin']);
    if nchans~=hdr.nChans && length(chanindx)==nchans
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      hdr.label = hdr.label(chanindx);
      hdr.nChans = length(chanindx);
    end
    % write the header file
    save(headerfile, 'hdr', '-v6');
    % write the data file
    [fid,message] = fopen(datafile,'wb','ieee-le');
    fwrite(fid, dat, 'double');
    fclose(fid);
    
  case 'fcdc_mysql'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write to a MySQL server listening somewhere else on the network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    db_open(filename);
    if ~isempty(hdr) && isempty(dat)
      % insert the header information into the database
      if db_blob
        % insert the structure into the database table as a binary blob
        db_insert_blob('fieldtrip.header', 'msg', hdr);
      else
        % make a structure with the same elements as the fields in the database table
        s             = struct;
        s.Fs          = hdr.Fs;           % sampling frequency
        s.nChans      = hdr.nChans;       % number of channels
        s.nSamples    = hdr.nSamples;     % number of samples per trial
        s.nSamplesPre = hdr.nSamplesPre;  % number of pre-trigger samples in each trial
        s.nTrials     = hdr.nTrials;      % number of trials
        s.label       = serialize(hdr.label); % FIXME this is not generic
        try
          s.msg = serialize(hdr); % FIXME this is not generic
        catch
          warning(lasterr);
        end
        db_insert('fieldtrip.header', s);
      end

    elseif isempty(hdr) && ~isempty(dat)
        dim = size(dat);
        if numel(dim)==2
          % ensure that the data dimensions correspond to ntrials X nchans X samples
          dim = [1 dim];
          dat = reshape(dat, dim);
        end
        ntrials = dim(1);
        for i=1:ntrials
          if db_blob
            % insert the data into the database table as a binary blob
            db_insert_blob('fieldtrip.data', 'msg', reshape(dat(i,:,:), dim(2:end)));
          else
            % create a structure with the same fields as the database table
            s = struct;
            s.nChans   = dim(2);
            s.nSamples = dim(3);
            try
              s.data = serialize(reshape(dat(i,:,:), dim(2:end)));
            catch
              warning(lasterr);
            end
            % insert the structure into the database
            db_insert('fieldtrip.data', s);
          end
        end
      
    else
      error('you should specify either the header or the data when writing to a MySQL database');
    end

  case 'plexon_nex'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single or mulitple channel Plexon NEX file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file, '.nex']);
    if nchans~=1
      error('only supported for single-channel data');
    end
    % construct a NEX structure with  the required parts of the header
    nex.hdr.VarHeader.Type       = 5; % continuous
    nex.hdr.VarHeader.Name       = hdr.label{1};
    nex.hdr.VarHeader.WFrequency = hdr.Fs;
    if isfield(hdr, 'FirstTimeStamp')
      nex.hdr.FileHeader.Frequency = hdr.Fs * hdr.TimeStampPerSample;
      nex.var.ts = hdr.FirstTimeStamp;
    else
      warning('no timestamp information available');
      nex.hdr.FileHeader.Frequency  = nan;
      nex.var.ts = nan;
    end
    nex.var.indx = 0;
    nex.var.dat  = dat;

    write_plexon_nex(filename, nex);

    if 0
      % the following code snippet can be used for testing
      nex2 = [];
      [nex2.var, nex2.hdr] = read_plexon_nex(filename, 'channel', 1);
    end

  case 'neuralynx_ncs'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single channel Neuralynx NCS file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [path, file, ext] = fileparts(filename);
    filename = fullfile(path, [file, '.ncs']);
    if nchans>1
      error('only supported for single-channel data');
    end

    if nchans~=hdr.nChans && length(chanindx)==nchans
      % assume that the header corresponds to the original multichannel
      % file and that the data represents a subset of channels
      % WARNING the AD channel index assumes that the data was read from a DMA or SDMA file
      % the first 17 channels contain status info, this number is zero-offset
      ADCHANNEL  = chanindx - 17 - 1;
      LABEL      = hdr.label{chanindx};
    elseif hdr.nChans==1
      ADCHANNEL  = -1;            % unknown
      LABEL      = hdr.label{1};  % single channel
    else
      error('cannot determine channel label');
    end

    FSAMPLE    = hdr.Fs;
    RECORDNSMP = 512;
    RECORDSIZE = 1044;

    % cut the downsampled LFP data into record-size pieces
    nrecords = ceil(nsamples/RECORDNSMP);
    fprintf('construct ncs with %d records\n', nrecords);

    % construct a ncs structure with all header details and the data in it
    ncs                = [];
    ncs.NumValidSamp   = ones(1,nrecords) * RECORDNSMP;   % except for the last block
    ncs.ChanNumber     = ones(1,nrecords) * ADCHANNEL;
    ncs.SampFreq       = ones(1,nrecords) * FSAMPLE;
    ncs.TimeStamp      = zeros(1,nrecords,'uint64');

    if rem(nsamples, RECORDNSMP)>0
      % the data length is not an integer number of records, pad the last record with zeros
      dat = cat(2, dat, zeros(nchans, nrecords*RECORDNSMP-nsamples));
      ncs.NumValidSamp(end) = rem(nsamples, RECORDNSMP);
    end

    ncs.dat = reshape(dat, RECORDNSMP, nrecords);

    for i=1:nrecords
      % timestamps should be 64 bit unsigned integers
      ncs.TimeStamp(i) = uint64(hdr.FirstTimeStamp) + uint64((i-1)*RECORDNSMP*hdr.TimeStampPerSample);
    end

    % add the elements that will go into the ascii header
    ncs.hdr.CheetahRev            = '4.23.0';
    ncs.hdr.NLX_Base_Class_Type   = 'CscAcqEnt';
    ncs.hdr.NLX_Base_Class_Name   = LABEL;
    ncs.hdr.RecordSize            = RECORDSIZE;
    ncs.hdr.ADChannel             = ADCHANNEL;
    ncs.hdr.SamplingFrequency     = FSAMPLE;

    % write it to a file
    fprintf('writing to %s\n', filename);
    write_neuralynx_ncs(filename, ncs);

    if 0
      % the following code snippet can be used for testing
      ncs2 = read_neuralynx_ncs(filename, 1, inf);
    end

  otherwise
    error('unsupported data format');
end % switch dataformat

