function [dat] = read_neuralynx_bin(filename, begsample, endsample);

% READ_NEURALYNX_BIN
%
% Use as
%   hdr = read_neuralynx_bin(filename)
% or
%   dat = read_neuralynx_bin(filename, begsample, endsample)
%
% This  is not a formal Neuralynx file format, but at the
% F.C. Donders Centre we use it in conjunction with Neuralynx,
% SPIKESPLITTING and SPIKEDOWNSAMPLE.

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: read_neuralynx_bin.m,v $
% Revision 1.2  2007/12/17 16:23:44  roboos
% fixed bug in jumping to correct begin sample
% added support for determining channel name from filename like this "dataset.chanlabel.bin"
% added support for old splitted dma files, which have an 8 byte header with the channel label and are always int32
%
% Revision 1.1  2007/12/12 16:28:42  roboos
% first implementation
%

needhdr = (nargin==1);
needdat = (nargin>=2);

% this is used for backward compatibility
oldformat = false;

% the first 8 bytes contain the header
fid    = fopen(filename, 'rb', 'ieee-le');
magic  = fread(fid, 8, 'char=>char')';

% the header describes the format of the subsequent samples
format = strtrim(magic(magic~=0));  % remove padding zeros and trailing spaces
if     strcmp(format, 'uint8')
  samplesize = 1;
elseif strcmp(format, 'int8')
  samplesize = 1;
elseif strcmp(format, 'uint16')
  samplesize = 2;
elseif strcmp(format, 'int16')
  samplesize = 2;
elseif strcmp(format, 'uint32')
  samplesize = 4;
elseif strcmp(format, 'int32')
  samplesize = 4;
elseif strcmp(format, 'uint64')
  samplesize = 8;
elseif strcmp(format, 'int64')
  samplesize = 8;
elseif strcmp(format, 'float32')
  samplesize = 4;
elseif strcmp(format, 'float64')
  samplesize = 8;
else
  warning('could not detect sample format, assuming old format file with ''int32''');
  oldformat  = true;
  samplesize = 4;
  format     = 'int32';
end

if oldformat
  % the first 8 bytes of the file contain the channel label (padded with spaces)
  label = strtrim(magic);
else
  % the filename is formatted like "dataset.chanlabel.bin"
  [p, f, x1] = fileparts(filename);
  [p, f, x2] = fileparts(f);
  if isempty(x2)
    warning('could not determine channel label');
    label = 'unknown';
  else
    label = x2(2:end);
  end
  clear p f x1 x2
end

% construct the header
hdr             = [];
hdr.Fs          = 32556;
hdr.nChans      = 1;
hdr.nSamples    = (filesize(filename)-8)/samplesize;
hdr.nSamplesPre = 0;
hdr.nTrials     = 1;
hdr.label       = {label};

if ~needdat
  % return only the header details
  dat = hdr;

else
  % read and return the data
  if begsample<1
    begsample = 1;
  end

  if isinf(endsample)
    endsample = hdr.nSamples;
  end

  fseek(fid, 8+(begsample-1)*samplesize, 'bof');   % skip to the beginning of the interesting data
  format = sprintf('%s=>%s', format, format);
  dat = fread(fid, [1 endsample-begsample+1], format);
  if length(dat)<(endsample-begsample+1)
    error('could not read the requested data');
  end
end % needdat

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the file size in bytes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [siz] = filesize(filename)
l = dir(filename);
if l.isdir
  error(sprintf('"%s" is not a file', filename));
end
siz = l.bytes;

