function [dat] = read_neuralynx_ds(dirname, hdr, begsample, endsample, chanindx);

% READ_NEURALYNX_DS reads multiple single-channel Neuralynx files that are
% all contained in a single directory. Each file is treated as a single
% channel of a combined multi-channel dataset.
%
% Use as
%   [hdr] = read_neuralynx_ds(dirname)
%   [dat] = read_neuralynx_ds(dirname, hdr, begsample, endsample, chanindx)
%
% A Neuralynx dataset consists of a directory containing separate files,
% one for each channel. All Neuralynx datafiles starts with a 16k header
% (in ascii format), followed by an arbitrary number of data records. The
% format of the data records depend on the type of data contained in the
% channel (e.g. continuous or spike data).
%
% To read the timestamps of spike waveforms (nse) or clustered spikes (nts),
% the header should contain the fields
%   hdr.FirstTimeStamp
%   hdr.TimeStampPerSample
% These can only be obtained from the corresponding simultaneous LFP
% and/or MUA recordings.
%
% See also READ_NEURALYNX_NCS, READ_NEURALYNX_NSE, READ_NEURALYNX_NTS

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: read_neuralynx_ds.m,v $
% Revision 1.3  2007/03/21 17:06:57  roboos
% updated the documentation
%
% Revision 1.2  2006/12/13 15:49:33  roboos
% many changes, cleaned up code, removed subfunctions, explicitely implement teh computation fo timestamp stuff using a mix of double and uint64 precision (using overloaded mex functions in @uint64 toolbox)
%
% Revision 1.1  2006/03/29 14:39:20  roboos
% New implementation, replaces the previous seperate read_neuralynx_header
% and read_neuralynx_data. Instead of reimplementing the reading code, this
% function calls the specialized subfunctions.
%

needhdr = (nargin==1);
needdat = (nargin>=2);

if needhdr

  % get the list of filenames
  ls = dir(dirname);
  ls = ls(~cell2mat({ls.isdir}));
  fname = {};
  for i=1:length(ls)
    fname{i} = fullfile(dirname, ls(i).name);
  end

  ftype = zeros(length(fname), 1);
  for i=1:length(fname)
    if     filetype(fname{i}, 'neuralynx_ncs')
      ftype(i) = 1;
    elseif filetype(fname{i}, 'neuralynx_nse')
      ftype(i) = 2;
    elseif filetype(fname{i}, 'neuralynx_nts')
      ftype(i) = 3;
    end
  end

  % only remember the filenames that are relevant
  fname = fname(ftype>0);
  ftype = ftype(ftype>0);

  if length(fname)==0
    error('the directory contains no supported files');
  end

  for i=1:length(fname)
    % this will only work if all files within a dataset return a similar header structure
    switch ftype(i)
      case 1
        orig(i) = read_neuralynx_ncs(fname{i});
      case 2
        orig(i) = read_neuralynx_nse(fname{i});
      case 3
        orig(i) = read_neuralynx_nts(fname{i});
      otherwise
        error('unsupported');
    end
  end

  % combine the information from the different files in a single header
  for i=1:length(orig)
    label{i}             = orig(i).hdr.NLX_Base_Class_Name;
    SamplingFrequency(i) = orig(i).hdr.SamplingFrequency;
    ADBitVolts(i)        = orig(i).hdr.ADBitVolts;
    NRecords(i)          = orig(i).NRecords;
    FirstTimeStamp(i)    = orig(i).TimeStamp(1);
    LastTimeStamp(i)     = orig(i).TimeStamp(end);
    % Note that the last timestamp corresponds with the first sample of the last
    % record and not with the last sample in the file.
  end

  for i=1:length(orig)
    % timestamps are measured in units of approximately 1us
    % in case of 32556 Hz sampling rate, there are approximately 30.7 timestamps per sample
    % in case of 1000 Hz sampling rate,  there are approximately 1000 timestamps per sample
    % note that the last timestamp in the original header corresponds with the
    % first sample of the last record, and not with the last sample
    switch ftype(i)
      case 1
        % ensure that the last timestamp refers to the last sample
        recordsize = 512; % each record contains 512 samples
      case 2
        recordsize = 32; % each record contains 32 samples
      case 3
        error('unsupported');
      case 4
        error('unsupported');
      otherwise
        error('unsupported');
    end
    % ensure that the last timestamp refers to the last sample
    TimeStampPerSample(i) = double(LastTimeStamp(i)-FirstTimeStamp(i))/((NRecords(i)-1)*recordsize);  % this should be in double precision, since it can be fractional
    LastTimeStamp(i)      = LastTimeStamp(i) + uint64((recordsize-1)*TimeStampPerSample(i));          % this should be in uint64 precision
  end % for length(orig)

  if any(SamplingFrequency~=SamplingFrequency(1))
    warning('not all channels have the same sampling rate');
  end

  if all(ftype==1) && any(FirstTimeStamp~=FirstTimeStamp(1))
    warning('not all channels start at the same time');
  end

  if all(ftype==1) && any(LastTimeStamp~=LastTimeStamp(1))
    warning('not all channels end at the same time');
  end

  if all(ftype==1) && any(NRecords~=NRecords(1))
    warning('the different channels do not have the same number of records');
  end

  % construct the header that applies to all channels combined
  hdr.nChans         = length(label);
  hdr.label          = label;
  hdr.filename       = fname;
  hdr.nTrials        = 1;                           % it is continuous
  hdr.Fs             = SamplingFrequency(1);
  hdr.nSamplesPre    = 0;                           % it is continuous

  if all(ftype==1)
    % these elements are only relevant for continously sampled channels
    hdr.nSamples           = NRecords(1) * 512;
    hdr.FirstTimeStamp     = FirstTimeStamp(1);
    hdr.LastTimeStamp      = LastTimeStamp(1);
    hdr.TimeStampPerSample = TimeStampPerSample(1);
  end

  % remember the original header details
  hdr.orig = orig(:);

  % return the header
  dat = hdr;

else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % read the data of the selected channels (i.e. files)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  if nargin<5
    % select all channels
    chanindx = 1:length(hdr.label);
  end
  nchan   = length(chanindx);
  nsample = endsample-begsample+1;
  dat     = zeros(nchan, nsample);

  for i=1:nchan
    thischan = chanindx(i);
    thisfile = hdr.filename{thischan};
    switch filetype(thisfile)
    case 'neuralynx_ncs'
      % determine the records that contain the sample numbers of the requested segment
      begrecord  = floor(begsample/512) + 1;
      endrecord  = ceil (endsample/512);
      ncs = read_neuralynx_ncs(thisfile, begrecord, endrecord);
      % copy the selected samples into the output
      begsel = begsample - (begrecord-1)*512;
      endsel = endsample - (begrecord-1)*512;
      dat(i,:) = ncs.dat(begsel:endsel);

    case 'neuralynx_nse'
      % read all spike waveforms and timestamps
      nse = read_neuralynx_nse(thisfile);
      % convert the timestamps into samples
      fprintf('%d timstamps\n', length(nse.TimeStamp));
      sample = double(nse.TimeStamp-hdr.FirstTimeStamp)/hdr.TimeStampPerSample + 1;
      sample = sample(sample>=begsample & sample<=endsample) - begsample + 1;
      dat(i,sample) = 1;

    case 'neuralynx_nts'
      % read all timestamps
      nts = read_neuralynx_nts(thisfile);
      % convert the timestamps into samples
      sample = double(nse.TimeStamp-hdr.FirstTimeStamp)/hdr.TimeStampPerSample + 1;
      sample = sample(sample>=begsample & sample<=endsample) - begsample + 1;
      dat(i,sample) = 1;

    end % switch filetype
  end % for nchan
end % reading data

