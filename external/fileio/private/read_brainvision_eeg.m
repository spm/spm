function [dat] = read_brainvision_eeg(filename, hdr, begsample, endsample);

% READ_BRAINVISION_EEG reads raw data from an EEG file
% and returns it as a Nchans x Nsamples matrix
%
% Use as
%   dat = read_brainvision_eeg(filename, hdr, begsample, endsample) 
% where the header should be first read using read_brainvision_vhdr
%
% See also READ_BRAINVISION_VHDR, READ_BRAINVISION_VMRK

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: read_brainvision_eeg.m,v $
% Revision 1.4  2007/06/13 08:08:19  roboos
% changed single & into &&
%
% Revision 1.3  2004/03/30 11:47:50  roberto
% dos->unix, fixed bug in multiplexed binary
%
% Revision 1.2  2004/03/30 08:22:19  roberto
% fixed bug due to renaming NumberOfChannels -> nChans
%
% Revision 1.1  2004/03/30 07:23:28  roberto
% added to CVS repository, copyrights added and implemented multiple file
% formats
%

if strcmpi(hdr.DataFormat, 'binary') && strcmpi(hdr.DataOrientation, 'multiplexed') && strcmpi(hdr.BinaryFormat, 'int_16')
  fid = fopen(filename, 'rb', 'ieee-le');
  fseek(fid, hdr.nChans*2*(begsample-1), 'cof');
  [dat, siz] = fread(fid, [hdr.nChans, (endsample-begsample+1)], 'int16');
  fclose(fid);
  % compute real microvolts using the calibration factor (resolution)
  res = sparse(diag(hdr.resolution));
  dat = res * dat;

elseif strcmpi(hdr.DataFormat, 'binary') && strcmpi(hdr.DataOrientation, 'multiplexed') && strcmpi(hdr.BinaryFormat, 'ieee_float_32')
  fid = fopen(filename, 'rb', 'ieee-le');
  fseek(fid, hdr.nChans*4*(begsample-1), 'cof');
  [dat, siz] = fread(fid, [hdr.nChans, (endsample-begsample+1)], 'float32');
  fclose(fid);

elseif strcmpi(hdr.DataFormat, 'binary') && strcmpi(hdr.DataOrientation, 'vectorized') && strcmpi(hdr.BinaryFormat, 'ieee_float_32')
  fid = fopen(filename, 'rb', 'ieee-le');
  fseek(fid, 0, 'eof');
  hdr.nSamples = ftell(fid)/(4*hdr.nChans);
  fseek(fid, 0, 'bof');
  numsamples = (endsample-begsample+1);
  for chan=1:hdr.nChans
    fseek(fid, (begsample-1)*4, 'cof');                 % skip the first N samples
    [tmp, siz] = fread(fid, numsamples, 'float32');     % read these samples
    fseek(fid, (hdr.nSamples-endsample)*4, 'cof');      % skip the last M samples
    dat(chan,:) = tmp(:)';
  end
  fclose(fid);

elseif strcmpi(hdr.DataFormat, 'ascii') && strcmpi(hdr.DataOrientation, 'multiplexed')
  fid = fopen(filename, 'rt');
  for line=1:(begsample-1)
    % read first lines and discard the data in them
    str = fgets(fid);
  end
  dat = zeros(endsample-begsample+1, hdr.nChans);
  for line=1:(endsample-begsample+1)
    str = fgets(fid);			% read a single line with Nchan samples
    str(find(str==',')) = '.';		% replace comma with point
    dat(line,:) = str2num(str);
  end
  fclose(fid);
  % transpose the data
  dat = dat';

elseif strcmpi(hdr.DataFormat, 'ascii') && strcmpi(hdr.DataOrientation, 'vectorized')
  % this is a very inefficient fileformat to read data from, since it requires to
  % read in all the samples of each channel and then select only the samples of interest
  fid = fopen(filename, 'rt');
  dat = zeros(hdr.nChans, endsample-begsample+1);
  for chan=1:hdr.nChans
    str = fgets(fid);			% read all samples of a single channel
    str(find(str==',')) = '.';		% replace comma with point
    tmp = str2num(str);
    dat(chan,:) = tmp(begsample:endsample);
  end
  fclose(fid);

else
  error('unsupported sub-fileformat');
end

