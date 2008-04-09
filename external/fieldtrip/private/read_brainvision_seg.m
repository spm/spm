function [dat] = read_brainvision_seg(filename, hdr, begsample, endsample);

% READ_BRAINVISION_SEG reads raw data from an segmented EEG file
% and returns it as a Nchans x Nsamples matrix.  The data is read
% as if it were continuous data, this function does not check whether
% the specified begin and endsample are withing the same segment.
%
% Use as
%   dat = read_brainvision_seg(filename, hdr, begsample, endsample) 
% where the header should be first read using read_brainvision_vhdr
%
% See also READ_BRAINVISION_VHDR, READ_BRAINVISION_VMRK

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: read_brainvision_seg.m,v $
% Revision 1.4  2008/04/09 10:11:23  roboos
% renamed nChans into the original form, converted from dos to unix
%
% Revision 1.3  2007/06/13 08:08:19  roboos
% changed single & into &&
%
% Revision 1.2  2004/03/30 08:22:19  roberto
% fixed bug due to renaming NumberOfChannels -> nChans
%
% Revision 1.1  2004/03/30 07:23:28  roberto
% added to CVS repository, copyrights added and implemented multiple file
% formats
%

if strcmpi(hdr.DataFormat, 'binary') && strcmpi(hdr.DataOrientation, 'multiplexed') && strcmpi(hdr.BinaryFormat, 'int_16')
  % this is a fileformat that I understand
  fid = fopen(filename, 'rb', 'ieee-le');
  fseek(fid, hdr.NumberOfChannels*2*(begsample-1), 'cof');
  [dat, siz] = fread(fid, [hdr.NumberOfChannels, (endsample-begsample+1)], 'int16');
  fclose(fid);
  % compute real microvolts using the calibration factor (resolution)
  res = sparse(diag(hdr.resolution));
  dat = res * dat;
else
  error('unsupported sub-fileformat');
end

