function [hdr] = read_header(filename, varargin)

% READ_HEADER is a wrapper around different EEG/MEG file importers
% directly supported formats are CTF, Neuromag, EEP, BrainVision,
% Neuroscan and Neuralynx.
%
% Use as
%   hdr = read_header(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'headerformat'   string
%
% This returns a header structure with the following elements
%   hdr.Fs           sampling frequency
%   hdr.nChans       number of channels
%   hdr.nSamples     number of samples per trial
%   hdr.nSamplesPre  number of pre-trigger samples in each trial
%   hdr.nTrials      number of trials
%   hdr.label        cell-array with labels of each channel
%
% For continuous data, nSamplesPre=0 and nTrials=1.
%
% Depending on the file format, additional header information can be
% returned in the hdr.orig subfield.
%
% See also READ_DATA, READ_EVENT, WRITE_DATA, WRITE_EVENT

% TODO channel renaming should be made a general option (see bham_bdf)

% Copyright (C) 2003-2007, Robert Oostenveld, F.C. Donders Centre
%
% $Log: read_header.m,v $
% Revision 1.37  2008/01/10 12:57:34  roboos
% give explicit errors with msgid FILEIO:Something
%
% Revision 1.36  2007/12/17 16:17:16  roboos
% updated some comments in the code, no functional change
%
% Revision 1.35  2007/12/17 08:24:28  roboos
% added support for nexstim_nxe, thanks to Vladimir
% the low-level code has not been tested by myself
%
% Revision 1.34  2007/12/12 16:50:15  roboos
% added support for neuralynx_bin
%
% Revision 1.33  2007/11/07 10:49:07  roboos
% cleaned up the reading and writing from/to mysql database, using db_xxx helper functions (see mysql directory)
%
% Revision 1.32  2007/11/05 17:01:21  roboos
% added implementation for fcdc_mysql
%
% Revision 1.31  2007/09/13 09:57:03  roboos
% use read_biosemi_bdf instead of openbdf/readbdf
% some small changes as sugegsted by the matlab editor (e.g. comma, semicolon)
%
% Revision 1.30  2007/08/01 12:24:40  roboos
% updated comments
%
% Revision 1.29  2007/08/01 09:57:01  roboos
% moved all code related to ctf shared memory to seperate functions
%
% Revision 1.28  2007/07/30 12:17:21  roboos
% ctf_shm: convert number of samples to double, otherwise problems with floor() in read_data
%
% Revision 1.27  2007/07/27 12:19:56  roboos
% added ctf_shm
%
% Revision 1.26  2007/07/19 14:49:30  roboos
% switched the default reader for nex files from read_nex_data to read_plexon_nex, the old one is still supported if explicitely mentioned as data/headerformat
%
% Revision 1.25  2007/07/04 13:20:51  roboos
% added support for egi_egis/egia, thanks to Joseph Dien
%
% Revision 1.24  2007/07/03 15:53:46  roboos
% switched from using Cristian Wienbruchs BTi toolbox to a new ascii header reading function (read_bti_m4d)
%
% Revision 1.23  2007/06/13 08:06:22  roboos
% updated help
%
% Revision 1.22  2007/06/06 12:39:43  roboos
% added try-catch for ctf2grad, to allow working with incomplete recordings on the new 275ch system
%
% Revision 1.21  2007/04/16 16:06:11  roboos
% add labels for neuroscan eeg format (thanks to Vladimir)
%
% Revision 1.20  2007/03/21 17:24:01  roboos
% added plexon_ds
%
% Revision 1.19  2007/03/19 17:07:19  roboos
% implemented an alternative reader for NEX (read_plexon_nex), the old implementation is still the default
%
% Revision 1.18  2007/02/21 09:54:21  roboos
% added timestamp details to header in case of neuralynx_ncs
%
% Revision 1.17  2007/01/10 17:29:54  roboos
% moved the fieldtrip specific handling of header information into read_header
%
% Revision 1.16  2007/01/09 09:38:40  roboos
% added spike channels for plexon_plx, moved plexon timestamp combination code to seperate function
%
% Revision 1.15  2007/01/04 17:13:19  roboos
% finished timestamp stuff for plexon_plx, only return the continuous channels that have data
%
% Revision 1.14  2007/01/04 12:20:26  roboos
% updated the neuralynx section to reflect the new reading functions and to use read_neuralynx_cds
% implemented recursive reading function for nested dataset directories
%
% Revision 1.13  2007/01/04 12:06:36  roboos
% added plexon_plx, not yet completely finished
%
% Revision 1.12  2006/12/04 10:37:10  roboos
% added support for ns_avg
% fixed bug in ns_eeg (hdr.nSamplesPre was negative)
%
% Revision 1.11  2006/10/09 15:39:51  roboos
% renamed BTi channels from 'MEGxxx' into 'Axxx'
%
% Revision 1.10  2006/09/18 21:48:33  roboos
% implemented support for fcdc_matbin, i.e. a dataset consisting of a matlab file with header and events and a seperate binary datafile
%
% Revision 1.9  2006/09/18 14:52:14  roboos
% implemented bti2grad and added it to header
%
% Revision 1.8  2006/09/18 14:22:54  roboos
% implemented support for 4D-BTi dataformat
%
% Revision 1.7  2006/09/13 11:01:01  roboos
% changed text in a error message
%
% Revision 1.6  2006/08/28 10:13:03  roboos
% use seperate filetype_check_extension function instead of subfunction, removed subfunction
%
% Revision 1.5  2006/06/22 15:07:39  roboos
% fiuxed bug in label assignment for tsl/tsh/ttl
%
% Revision 1.4  2006/06/22 13:51:47  roboos
% fixed typo in code: datatype instead of datatyppe, thanks to Thilo
%
% Revision 1.3  2006/06/22 07:53:46  roboos
% remember the original neuroscan cnt header
%
% Revision 1.2  2006/06/19 10:32:16  roboos
% added documentation
%
% Revision 1.1  2006/06/07 09:32:20  roboos
% new implementation based on the read_fcdc_xxx functions, now with
% variable (key-val) input arguments, changed the control structure
% in the rpobram (switch instead of ifs), allow the user to specify
% the file format, allow the user to specify either a sample selection
% or a block selection. The reading functionality should not have
% changed compared to the read_fcdc_xxx versions.
%

persistent db_blob       % for fcdc_mysql
if isempty(db_blob)
  db_blob = 0;
end

% test whether the file or directory exists
if ~exist(filename, 'file') && ~strcmp(filetype(filename), 'ctf_shm') && ~strcmp(filetype(filename), 'fcdc_mysql')
  error('FILEIO:InvalidFileName', 'file or directory ''%s'' does not exist', filename);
end

% get the options
headerformat = keyval('headerformat',  varargin);

% determine the filetype
if isempty(headerformat)
  headerformat = filetype(filename);
end

% start with an empty header
hdr = [];

switch headerformat
  case '4d_pdf'
    datafile   = filename;
    headerfile = [datafile '.m4d'];
    sensorfile = [datafile '.xyz'];
  case {'4d_m4d', '4d_xyz'}
    datafile   = filename(1:(end-4)); % remove the extension
    headerfile = [datafile '.m4d'];
    sensorfile = [datafile '.xyz'];
  case 'ctf_ds'
    % convert CTF filename into filenames
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(filename, [file '.res4']);
    datafile   = fullfile(filename, [file '.meg4']);
  case 'ctf_meg4'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.res4']);
    datafile   = fullfile(path, [file '.meg4']);
  case 'ctf_res4'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.res4']);
    datafile   = fullfile(path, [file '.meg4']);
  case 'brainvision_vhdr'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.vhdr']);
    if exist(fullfile(path, [file '.eeg']))
      datafile   = fullfile(path, [file '.eeg']);
    elseif exist(fullfile(path, [file '.seg']))
      datafile   = fullfile(path, [file '.seg']);
    elseif exist(fullfile(path, [file '.dat']))
      datafile   = fullfile(path, [file '.dat']);
    end
  case 'brainvision_eeg'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.vhdr']);
    datafile   = fullfile(path, [file '.eeg']);
  case 'brainvision_seg'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.vhdr']);
    datafile   = fullfile(path, [file '.seg']);
  case 'brainvision_dat'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.vhdr']);
    datafile   = fullfile(path, [file '.dat']);
  case 'fcdc_matbin'
    [path, file, ext] = fileparts(filename);
    headerfile = fullfile(path, [file '.mat']);
    datafile   = fullfile(path, [file '.bin']);
  otherwise
    % convert filename into filenames, assume that the header and data are the same
    datafile   = filename;
    headerfile = filename;
end

if ~strcmp(filename, headerfile)
  filename     = headerfile;                % this function will read the header
  headerformat = filetype(filename);        % update the filetype
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the data with the low-level reading function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch headerformat

  case {'4d_pdf', '4d_m4d', '4d_xyz'}
    orig = read_bti_m4d(filename);
    hdr.Fs          = orig.SampleFrequency;
    hdr.nChans      = orig.TotalChannels;
    hdr.nSamples    = orig.SlicesPerEpoch;
    hdr.nSamplesPre = round(-orig.FirstLatency*orig.SampleFrequency);
    hdr.nTrials     = orig.TotalEpochs;
    hdr.label       = orig.ChannelOrder(:);
    hdr.grad        = orig.grad;
    % remember original header details
    hdr.orig = orig;

  case 'besa_avr'
    orig = read_besa_avr(filename);
    hdr.Fs          = 1000/orig.di;
    hdr.nChans      = size(orig.data,1);
    hdr.nSamples    = size(orig.data,2);
    hdr.nSamplesPre = -(hdr.Fs * orig.tsb/1000);   % convert from ms to samples
    hdr.nTrials     = 1;
    if isfield(orig, 'label') && iscell(orig.label)
      hdr.label = orig.label;
    elseif isfield(orig, 'label') && ischar(orig.label)
      hdr.label = tokenize(orig.label, ' ');
    else
      for i=1:hdr.nChans
        warning('creating fake channel names for besa_avr');
        hdr.label{i} = sprintf('%03d', i);
      end
    end

  case 'besa_swf'
    orig = read_besa_swf(filename);
    hdr.Fs          = 1000/orig.di;
    hdr.nChans      = size(orig.data,1);
    hdr.nSamples    = size(orig.data,2);
    hdr.nSamplesPre = -(hdr.Fs * orig.tsb/1000);   % convert from ms to samples
    hdr.nTrials     = 1;
    hdr.label       = orig.label;
    
  case {'biosemi_bdf', 'bham_bdf'}
    hdr = read_biosemi_bdf(filename);
    if any(diff(hdr.orig.SampleRate))
      error('channels with different sampling rate not supported');
    end
    if filetype(filename, 'bham_bdf')
      % TODO channel renaming should be made a general option
      % this is for the Biosemi system used at the University of Birmingham
      labelold = { 'A1' 'A2' 'A3' 'A4' 'A5' 'A6' 'A7' 'A8' 'A9' 'A10' 'A11' 'A12' 'A13' 'A14' 'A15' 'A16' 'A17' 'A18' 'A19' 'A20' 'A21' 'A22' 'A23' 'A24' 'A25' 'A26' 'A27' 'A28' 'A29' 'A30' 'A31' 'A32' 'B1' 'B2' 'B3' 'B4' 'B5' 'B6' 'B7' 'B8' 'B9' 'B10' 'B11' 'B12' 'B13' 'B14' 'B15' 'B16' 'B17' 'B18' 'B19' 'B20' 'B21' 'B22' 'B23' 'B24' 'B25' 'B26' 'B27' 'B28' 'B29' 'B30' 'B31' 'B32' 'C1' 'C2' 'C3' 'C4' 'C5' 'C6' 'C7' 'C8' 'C9' 'C10' 'C11' 'C12' 'C13' 'C14' 'C15' 'C16' 'C17' 'C18' 'C19' 'C20' 'C21' 'C22' 'C23' 'C24' 'C25' 'C26' 'C27' 'C28' 'C29' 'C30' 'C31' 'C32' 'D1' 'D2' 'D3' 'D4' 'D5' 'D6' 'D7' 'D8' 'D9' 'D10' 'D11' 'D12' 'D13' 'D14' 'D15' 'D16' 'D17' 'D18' 'D19' 'D20' 'D21' 'D22' 'D23' 'D24' 'D25' 'D26' 'D27' 'D28' 'D29' 'D30' 'D31' 'D32' 'EXG1' 'EXG2' 'EXG3' 'EXG4' 'EXG5' 'EXG6' 'EXG7' 'EXG8' 'Status'};
      labelnew = { 'P9' 'PPO9h' 'PO7' 'PPO5h' 'PPO3h' 'PO5h' 'POO9h' 'PO9' 'I1' 'OI1h' 'O1' 'POO1' 'PO3h' 'PPO1h' 'PPO2h' 'POz' 'Oz' 'Iz' 'I2' 'OI2h' 'O2' 'POO2' 'PO4h' 'PPO4h' 'PO6h' 'POO10h' 'PO10' 'PO8' 'PPO6h' 'PPO10h' 'P10' 'P8' 'TPP9h' 'TP7' 'TTP7h' 'CP5' 'TPP7h' 'P7' 'P5' 'CPP5h' 'CCP5h' 'CP3' 'P3' 'CPP3h' 'CCP3h' 'CP1' 'P1' 'Pz' 'CPP1h' 'CPz' 'CPP2h' 'P2' 'CPP4h' 'CP2' 'CCP4h' 'CP4' 'P4' 'P6' 'CPP6h' 'CCP6h' 'CP6' 'TPP8h' 'TP8' 'TPP10h' 'T7' 'FTT7h' 'FT7' 'FC5' 'FCC5h' 'C5' 'C3' 'FCC3h' 'FC3' 'FC1' 'C1' 'CCP1h' 'Cz' 'FCC1h' 'FCz' 'FFC1h' 'Fz' 'FFC2h' 'FC2' 'FCC2h' 'CCP2h' 'C2' 'C4' 'FCC4h' 'FC4' 'FC6' 'FCC6h' 'C6' 'TTP8h' 'T8' 'FTT8h' 'FT8' 'FT9' 'FFT9h' 'F7' 'FFT7h' 'FFC5h' 'F5' 'AFF7h' 'AF7' 'AF5h' 'AFF5h' 'F3' 'FFC3h' 'F1' 'AF3h' 'Fp1' 'Fpz' 'Fp2' 'AFz' 'AF4h' 'F2' 'FFC4h' 'F4' 'AFF6h' 'AF6h' 'AF8' 'AFF8h' 'F6' 'FFC6h' 'FFT8h' 'F8' 'FFT10h' 'FT10'};
      % rename the channel labels
      for i=1:length(labelnew)
        chan = strmatch(labelold(i), hdr.label, 'exact');
        hdr.label(chan) = labelnew(chan);
      end
    end

  case {'biosemi_old'}
    % this uses the openbdf and readbdf functions that I copied from the EEGLAB toolbox
    orig = openbdf(filename);
    if any(orig.Head.SampleRate~=orig.Head.SampleRate(1))
      error('channels with different sampling rate not supported');
    end
    hdr.Fs          = orig.Head.SampleRate(1);
    hdr.nChans      = orig.Head.NS;
    hdr.label       = cellstr(orig.Head.Label);
    % it is continuous data, therefore append all records in one trial
    hdr.nSamples    = orig.Head.NRec * orig.Head.Dur * orig.Head.SampleRate(1);
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1;
    hdr.orig        = orig;
    % close the file between seperate read operations
    fclose(orig.Head.FILE.FID);

  case {'brainvision_vhdr', 'brainvision_seg', 'brainvision_eeg', 'brainvision_dat'}
    hdr = read_brainvision_vhdr(filename);

  case 'ced_son'
    % check that the required low-level toolbox is available
    hastoolbox('neuroshare', 1);
    % use the reading function supplied by Gijs van Elswijk
    orig = read_ced_son(filename,'readevents','no','readdata','no');
    orig = orig.header;
    % In Spike2, channels can have different sampling rates, units, length
    % etc. etc. Here, channels need to have to same properties.
    if length(unique([orig.samplerate]))>1,
      error('channels with different sampling rates are not supported');
    else
      hdr.Fs   = orig(1).samplerate;
    end;
    hdr.nChans = length(orig);
    % nsamples of the channel with least samples
    hdr.nSamples    = min([orig.nsamples]);
    hdr.nSamplesPre = 0;
    % only continuous data supported
    if sum(strcmpi({orig.mode},'continuous')) < hdr.nChans,
      error('not all channels contain continuous data');
    else,
      hdr.nTrials = 1;
    end;
    hdr.label = {orig.label};

  case {'ctf_ds', 'ctf_meg4', 'ctf_res4', 'read_ctf_res4'} % the default reader for CTF is read_ctf_res4
    if strcmp(headerformat, 'ctf_ds')
      [p, f, x] = fileparts(filename);
      filename = fullfile(filename, [f '.res4']);
    elseif strcmp(headerformat, 'ctf_meg4')
      [p, f, x] = fileparts(filename);
      filename = fullfile(p, [f '.res4']);
    end
    % read it using the default CTF importer that originates from CTF and the FCDC
    orig = read_ctf_res4(filename);
    hdr.Fs           = orig.Fs;
    hdr.nChans       = orig.nChans;
    hdr.nSamples     = orig.nSamples;
    hdr.nSamplesPre  = orig.nSamplesPre;
    hdr.nTrials      = orig.nTrials;
    hdr.label        = orig.label;
    % add a gradiometer structure for forward and inverse modelling
    try
      hdr.grad = ctf2grad(orig);
    catch
      warning('could not construct gradiometer definition from the header');
    end
    % add the original header details
    hdr.orig = orig;

  case 'ctf_read_res4'
    % check that the required low-level toolbos ix available
    hastoolbox('eegsf', 1);
    if strcmp(headerformat, 'ctf_res4')
      [p, f, x] = fileparts(filename);
      filename = p;
    elseif strcmp(headerformat, 'ctf_meg4')
      [p, f, x] = fileparts(filename);
      filename = p;
    end
    % read it using the CTF importer from the NIH and Daren Weber
    orig = ctf_read_res4(filename, 0);
    % convert the header into a structure that FieldTrip understands
    hdr              = [];
    hdr.Fs           = orig.setup.sample_rate;
    hdr.nChans       = length(orig.sensor.info);
    hdr.nSamples     = orig.setup.number_samples;
    hdr.nSamplesPre  = orig.setup.pretrigger_samples;
    hdr.nTrials      = orig.setup.number_trials;
    for i=1:length(orig.sensor.info)
      hdr.label{i}   = orig.sensor.info(i).label;
    end
    hdr.label        = hdr.label(:);
    % add a gradiometer structure for forward and inverse modelling
    hdr.grad         = nimh2grad(orig);
    % add the original header details
    hdr.orig = orig;

  case 'ctf_shm'
    % contact Robert Oostenveld if you are interested in real-time acquisition on the CTF system
    % read the header information from shared memory
    hdr = read_shm_header(filename);

  case 'eep_avr'
    % check that the required low-level toolbox is available
    hastoolbox('eeprobe', 1);
    % read the whole average and keep only header info (it is a bit silly, but the easiest to do here)
    hdr = read_eep_avr(filename);
    hdr.Fs          = hdr.rate;
    hdr.nChans      = size(hdr.data,1);
    hdr.nSamples    = size(hdr.data,2);
    hdr.nSamplesPre = hdr.xmin*hdr.rate/1000;
    hdr.nTrials     = 1;		% it can always be interpreted as continuous data
    % remove the data and variance if present
    hdr = rmfield(hdr, 'data');
    try, hdr = rmfield(hdr, 'variance'); end

  case 'eep_cnt'
    % check that the required low-level toolbox is available
    hastoolbox('eeprobe', 1);
    % read the first sample from the continous data, which will also return the header
    hdr = read_eep_cnt(filename, 1, 1);
    hdr.Fs          = hdr.rate;
    hdr.nSamples    = hdr.nsample;
    hdr.nSamplesPre = 0;
    hdr.nChans      = hdr.nchan;
    hdr.nTrials     = 1;		% it can always be interpreted as continuous data

  case 'egi_egia'
    [fhdr,chdr,ename,cnames,fcom,ftext] = read_egis_header(filename);
    [p, f, x]       = fileparts(filename);

    if any(chdr(:,4)-chdr(1,4))
      error('Sample rate not the same for all cells.');
    end;

    hdr.Fs          = chdr(1,4); %making assumption that sample rate is same for all cells
    hdr.nChans      = fhdr(19);
    for i = 1:hdr.nChans
      hdr.label{i}  = ['e' num2str(i)];
    end;
    hdr.nTrials     = fhdr(11)*fhdr(18); %number of trials is numSubjects * numCells
    hdr.nSamplesPre = ceil(fhdr(14)/(1000/hdr.Fs));

    if any(chdr(:,3)-chdr(1,3))
      error('Number of samples not the same for all cells.');
    end;

    hdr.nSamples    = chdr(1,3); %making assumption that number of samples is same for all cells

    % remember the original header details
    hdr.orig.fhdr   = fhdr;
    hdr.orig.chdr   = chdr;
    hdr.orig.ename  = ename;
    hdr.orig.cnames = cnames;
    hdr.orig.fcom   = fcom;
    hdr.orig.ftext  = ftext;

  case 'egi_egis'
    [fhdr,chdr,ename,cnames,fcom,ftext] = read_egis_header(filename);
    [p, f, x]       = fileparts(filename);

    if any(chdr(:,4)-chdr(1,4))
      error('Sample rate not the same for all cells.');
    end;

    hdr.Fs          = chdr(1,4); %making assumption that sample rate is same for all cells
    hdr.nChans      = fhdr(19);
    for i = 1:hdr.nChans
      hdr.label{i}  = ['e' num2str(i)];
    end;
    hdr.nTrials     = sum(chdr(:,2));
    hdr.nSamplesPre = ceil(fhdr(14)/(1000/hdr.Fs));  
    % assuming that a utility was used to insert the correct baseline
    % duration into the header since it is normally absent. This slot is
    % actually allocated to the age of the subject, although NetStation
    % does not use it when generating an EGIS session file.

    if hdr.nSamplesPre == 0
      hdr.nSamplesPre = 1; %If baseline was left as zero, then change to "1" to avoid possible issues with software expecting a non-zero baseline.
    end;

    if any(chdr(:,3)-chdr(1,3))
      error('Number of samples not the same for all cells.');
    end;

    hdr.nSamples    = chdr(1,3); %making assumption that number of samples is same for all cells

    % remember the original header details
    hdr.orig.fhdr   = fhdr;
    hdr.orig.chdr   = chdr;
    hdr.orig.ename  = ename;
    hdr.orig.cnames = cnames;
    hdr.orig.fcom   = fcom;
    hdr.orig.ftext  = ftext;

  case 'fcdc_matbin'
    % this is multiplexed data in a *.bin file, accompanied by a matlab file containing the header
    load(headerfile, 'hdr');

  case 'fcdc_mysql'
    % read from a MySQL server listening somewhere else on the network
    db_open(filename);
    if db_blob
      hdr = db_select_blob('fieldtrip.header', 'msg', 1);
    else
      hdr = db_select('fieldtrip.header', {'nChans', 'nSamples', 'nSamplesPre', 'Fs', 'label'}, 1);
      hdr.label = eval(hdr.label);  % FIXME this is not generic
    end

  case {'mpi_ds', 'mpi_dap'}
    hdr = read_mpi_ds(filename);

  case 'neuralynx_dma'
    hdr = read_neuralynx_dma(filename);

  case 'neuralynx_sdma'
    hdr = read_neuralynx_sdma(filename);

  case 'neuralynx_ncs'
    ncs = read_neuralynx_ncs(filename, 1, 0);
    [p, f, x]       = fileparts(filename);
    hdr.Fs          = ncs.hdr.SamplingFrequency;
    hdr.label       = {f};
    hdr.nChans      = 1;
    hdr.nTrials     = 1;
    hdr.nSamplesPre = 0;
    hdr.nSamples    = ncs.NRecords * 512;
    hdr.orig        = ncs.hdr;
    FirstTimeStamp  = ncs.hdr.FirstTimeStamp;  % this is the first timestamp of the first block
    LastTimeStamp   = ncs.hdr.LastTimeStamp;   % this is the first timestamp of the last block, i.e. not the timestamp of the last sample
    hdr.TimeStampPerSample = double(LastTimeStamp - FirstTimeStamp) ./ ((ncs.NRecords-1)*512);
    hdr.FirstTimeStamp     = FirstTimeStamp;
    hdr.LastTimeStamp      = LastTimeStamp + uint64((512-1)*hdr.TimeStampPerSample);

  case 'neuralynx_nse'
    nse = read_neuralynx_nse(filename, 1, 0);
    [p, f, x]       = fileparts(filename);
    hdr.Fs          = nse.hdr.SamplingFrequency;
    hdr.label       = {f};
    hdr.nChans      = 1;
    hdr.nTrials     = nse.NRecords;  % each record contains one waveform
    hdr.nSamples    = 32;            % there are 32 samples in each waveform
    hdr.nSamplesPre = 0;
    hdr.orig        = nse.hdr;
    % FIXME add hdr.FirstTimeStamp and hdr.LastTimeStamp

  case {'neuralynx_ttl', 'neuralynx_tsl', 'neuralynx_tsh'}
    % these are hardcoded, they contain an 8-byte header and int32 values for a single channel 
    % FIXME this should be done similar as neuralynx_bin, i.e. move the hdr into the function
    hdr             = [];
    hdr.Fs          = 32556;
    hdr.nChans      = 1;
    hdr.nSamples    = (filesize(filename)-8)/4;
    hdr.nSamplesPre = 1;
    hdr.nTrials     = 1;
    hdr.label       = {headerformat((end-3):end)};

  case 'neuralynx_bin'
    hdr = read_neuralynx_bin(filename);
    
  case 'neuralynx_ds'
    hdr = read_neuralynx_ds(filename);

  case 'neuralynx_cds'
    hdr = read_neuralynx_cds(filename);

  case 'nexstim_nxe'
    hdr = read_nexstim_nxe(filename);

  case 'neuromag_fif'
    % check that the required low-level toolbox is available
    hastoolbox('meg-pd', 1);
    % add a gradiometer structure for forward and inverse modelling
    hdr.grad = fif2grad(filename);
    [hdr.label, type, number] = channames(filename);
    rawdata('any',filename);
    hdr.Fs = rawdata('sf');
    hdr.nTrials = 0;
    rawdata('goto', 0);
    [buf, status] = rawdata('next');
    hdr.nChans   = size(buf,1);
    hdr.nSamples = size(buf,2);
    while strcmp(status, 'ok')
      hdr.nTrials = hdr.nTrials + 1;
      [buf, status] = rawdata('next');
    end
    % I don't know how to get this out of the file
    hdr.nSamplesPre  = 0;
    % this would give some information that I don't know how to use
    % [range, cal] = rawdata('range');
    % this would give the current latency
    % rawdata('t');
    rawdata('close');

  case 'ns_avg'
    orig = read_ns_hdr(filename);
    % do some reformatting/renaming of the header items
    hdr.Fs          = orig.rate;
    hdr.nSamples    = orig.npnt;
    hdr.nSamplesPre = round(-orig.rate*orig.xmin/1000);
    hdr.nChans      = orig.nchan;
    hdr.nTrials     = 1; % the number of trials in this datafile is only one, i.e. the average
    % remember the original header details
    hdr.orig = orig;

  case 'ns_cnt'
    % read_ns_cnt originates from the EEGLAB package (loadcnt.m) but is
    % an old version since the new version is not compatible any more
    orig = read_ns_cnt(filename, 'ldheaderonly', 1);
    % do some reformatting/renaming of the header items
    hdr.Fs          = orig.rate;
    hdr.nChans      = orig.nchannels;
    hdr.nSamples    = orig.nsamples;
    hdr.nSamplesPre = 0;
    hdr.nTrials     = 1;
    for i=1:hdr.nChans
      hdr.label{i} = deblank(orig.chan.names(i,:));
    end
    % remember the original header details
    hdr.orig = orig;

  case 'ns_eeg'
    orig = read_ns_hdr(filename);
    % do some reformatting/renaming of the header items
    hdr.label       = orig.label;
    hdr.Fs          = orig.rate;
    hdr.nSamples    = orig.npnt;
    hdr.nSamplesPre = round(-orig.rate*orig.xmin/1000);
    hdr.nChans      = orig.nchan;
    hdr.nTrials     = orig.nsweeps;
    % remember the original header details
    hdr.orig = orig;

  case 'plexon_ds'
    hdr = read_plexon_ds(filename);

  case 'plexon_ddt'
    orig = read_plexon_ddt(filename);
    hdr.nChans      = orig.NChannels;
    hdr.Fs          = orig.Freq;
    hdr.nSamples    = orig.NSamples;
    hdr.nSamplesPre = 0;      % continuous
    hdr.nTrials     = 1;      % continuous
    for i=1:hdr.nChans
      % create fake labels
      hdr.label{i} = sprintf('%03d', i);
    end
    % also remember the original header
    hdr.orig        = orig;

  case {'read_nex_data'} % this is an alternative reader for nex files
    orig = read_nex_header(filename);
    % assign the obligatory items to the output FCDC header
    numsmp = cell2mat({orig.varheader.numsmp});
    adindx = find(cell2mat({orig.varheader.typ})==5);
    if isempty(adindx)
      error('file does not contain continuous channels');
    end
    hdr.nChans      = length(orig.varheader);
    hdr.Fs          = orig.varheader(adindx(1)).wfrequency;     % take the sampling frequency from the first A/D channel
    hdr.nSamples    = max(numsmp(adindx));                      % take the number of samples from the longest A/D channel
    hdr.nTrials     = 1;                                        % it can always be interpreted as continuous data
    hdr.nSamplesPre = 0;                                        % and therefore it is not trial based
    for i=1:hdr.nChans
      hdr.label{i} = deblank(char(orig.varheader(i).nam));
    end
    hdr.label = hdr.label(:);
    % also remember the original header details
    hdr.orig = orig;

  case {'read_plexon_nex' 'plexon_nex'} % this is the default reader for nex files
    orig = read_plexon_nex(filename);
    numsmp = cell2mat({orig.VarHeader.NPointsWave});
    adindx = find(cell2mat({orig.VarHeader.Type})==5);
    if isempty(adindx)
      error('file does not contain continuous channels');
    end
    hdr.nChans      = length(orig.VarHeader);
    hdr.Fs          = orig.VarHeader(adindx(1)).WFrequency;     % take the sampling frequency from the first A/D channel
    hdr.nSamples    = max(numsmp(adindx));                      % take the number of samples from the longest A/D channel
    hdr.nTrials     = 1;                                        % it can always be interpreted as continuous data
    hdr.nSamplesPre = 0;                                        % and therefore it is not trial based
    for i=1:hdr.nChans
      hdr.label{i} = deblank(char(orig.VarHeader(i).Name));
    end
    hdr.label = hdr.label(:);
    hdr.FirstTimeStamp     = orig.FileHeader.Beg;
    hdr.TimeStampPerSample = orig.FileHeader.Frequency ./ hdr.Fs;
    % also remember the original header details
    hdr.orig = orig;

  case 'plexon_plx'
    orig = read_plexon_plx(filename);
    if orig.NumSlowChannels==0
      error('file does not contain continuous channels');
    end
    fsample = [orig.SlowChannelHeader.ADFreq];
    if any(fsample~=fsample(1))
      error('different sampling rates in continuous data not supported');
    end
    for i=1:length(orig.SlowChannelHeader)
      label{i} = deblank(orig.SlowChannelHeader(i).Name);
    end
    % continuous channels don't always contain data, remove the empty ones
    sel  = [orig.DataBlockHeader.Type]==5;  % continuous
    chan = [orig.DataBlockHeader.Channel];
    for i=1:length(label)
      chansel(i) = any(chan(sel)==orig.SlowChannelHeader(i).Channel);
    end
    chansel = find(chansel); % this is required for timestamp selection
    label = label(chansel);
    % only the continuous channels are returned as visible
    hdr.nChans      = length(label);
    hdr.Fs          = fsample(1);
    hdr.label       = label;
    % also remember the original header
    hdr.orig        = orig;

    % select the first continuous channel that has data
    sel = ([orig.DataBlockHeader.Type]==5 & [orig.DataBlockHeader.Channel]==orig.SlowChannelHeader(chansel(1)).Channel);
    % get the timestamps that correspond with the continuous data
    tsl = [orig.DataBlockHeader(sel).TimeStamp]';
    tsh = [orig.DataBlockHeader(sel).UpperByteOf5ByteTimestamp]';
    ts  = timestamp_plexon(tsl, tsh);  % use helper function, this returns an uint64 array

    % determine the number of samples in the continuous channels
    num = [orig.DataBlockHeader(sel).NumberOfWordsInWaveform];
    hdr.nSamples    = sum(num);
    hdr.nSamplesPre = 0;      % continuous
    hdr.nTrials     = 1;      % continuous

    % the timestamps indicate the beginning of each block, hence the timestamp of the last block corresponds with the end of the previous block
    hdr.TimeStampPerSample = double(ts(end)-ts(1))/sum(num(1:(end-1)));
    hdr.FirstTimeStamp     = ts(1);                                                %  the timestamp of the first continuous sample
    hdr.LastTimeStamp      = ts(end) + uint64(hdr.TimeStampPerSample * num(end));  % the timestamp of the last sample (not the beginning of the last block)

    % also make the spike channels visible
    for i=1:length(orig.ChannelHeader)
      hdr.label{end+1} = deblank(orig.ChannelHeader(i).Name);
    end
    hdr.label = hdr.label(:);
    hdr.nChans = length(hdr.label);

  case {'yokogawa_ave', 'yokogawa_con', 'yokogawa_raw'}
    % chek that the required low-level toolbox is available
    hastoolbox('yokogawa', 1);
    hdr = read_yokogawa_header(filename);
    % add a gradiometer structure for forward and inverse modelling
    hdr.grad = yokogawa2grad(hdr);

  otherwise
    error('unsupported header format');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the file size in bytes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [siz] = filesize(filename)
l = dir(filename);
if l.isdir
  error(sprintf('"%s" is not a file', filename));
end
siz = l.bytes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the file size in bytes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hdr] = recursive_read_header(filename)
[p, f, x] = fileparts(filename);
ls = dir(filename);
ls = ls(~strcmp({ls.name}, '.'));  % exclude this directory
ls = ls(~strcmp({ls.name}, '..')); % exclude parent directory
for i=1:length(ls)
  % make sure that the directory listing includes the complete path
  ls(i).name = fullfile(filename, ls(i).name);
end
lst = {ls.name};
hdr = cell(size(lst));
sel = zeros(size(lst));
for i=1:length(lst)
  % read the header of each individual file
  try
    thishdr = read_header(lst{i});
    if isstruct(thishdr)
      thishdr.filename = lst{i};
    end
  catch
    thishdr = [];
    warning(lasterr);
    fprintf('while reading %s\n\n', lst{i});
  end
  if ~isempty(thishdr)
    hdr{i} = thishdr;
    sel(i) = true;
  else
    sel(i) = false;
  end
end
sel = logical(sel(:));
lst = lst(sel);
hdr = hdr(sel);
tmp = {};
for i=1:length(hdr)
  if isstruct(hdr{i})
    tmp = cat(1, tmp, hdr(i));
  elseif iscell(hdr{i})
    tmp = cat(1, tmp, hdr{i}{:});
  end
end
hdr = tmp;
