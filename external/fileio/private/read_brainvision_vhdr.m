function [hdr] = read_brainvision_vhdr(filename);

% READ_BRAINVISION_VHDR reads the known items from the BrainVision EEG
% header file and returns them in a structure
%
% Use as
%   hdr = read_brainvision_vhdr(filename) 
%
% See also READ_BRAINVISION_EEG, READ_BRAINVISION_VMRK

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: read_brainvision_vhdr.m,v $
% Revision 1.2  2004/03/30 11:23:50  roberto
% fixed bug in cutting channel info into pieces, added local tokenize function
%
% Revision 1.1  2004/03/30 07:23:28  roberto
% added to CVS repository, copyrights added and implemented multiple file
% formats
%

hdr.DataFile         = read_asa(filename, 'DataFile=', '%s');
hdr.MarkerFile       = read_asa(filename, 'MarkerFile=', '%s');
hdr.DataFormat       = read_asa(filename, 'DataFormat=', '%s');
hdr.DataOrientation  = read_asa(filename, 'DataOrientation=', '%s');
hdr.BinaryFormat     = read_asa(filename, 'BinaryFormat=', '%s');
hdr.nChans           = read_asa(filename, 'NumberOfChannels=', '%d');
hdr.SamplingInterval = read_asa(filename, 'SamplingInterval=', '%f');	% microseconds

if ~isempty(hdr.nChans)
  for i=1:hdr.nChans
    chan_str  = sprintf('Ch%d=', i);
    chan_info = read_asa(filename, chan_str, '%s');
    t = tokenize(chan_info, ',');
    hdr.label{i} = t{1};
    hdr.reference{i} = t{2};
    resolution = str2num(t{3});          % in microvolt
    if ~isempty(resolution)
      hdr.resolution(i) = resolution;
    else
      hdr.resolution(i) = nan;
    end
  end
end

% compute the sampling rate in Hz
hdr.Fs = 1e6/(hdr.SamplingInterval);

% the number of samples is unkown
hdr.nSamples = Inf;

% the number of trials is unkown, assume continuous data
hdr.nTrials = 1;

% ensure that the labels are in a column
hdr.label      = hdr.label(:);
hdr.reference  = hdr.reference(:);
hdr.resolution = hdr.resolution(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction TOKENIZE cuts a string into pieces, returning a cell array
function t = tokenize(str, sep)
t = {};
f = find(str==sep);
if ~isempty(f)
  t{end+1} = str(1:(f(1)-1));
  for i=1:(length(f)-1)
    t{end+1} = str((f(i)+1):(f(i+1)-1));
  end
  t{end+1} = str((f(end)+1):end);
end
