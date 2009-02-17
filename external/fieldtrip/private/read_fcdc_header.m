function [hdr] = read_fcdc_header(filename)

% READ_FCDC_HEADER is a wrapper around different EEG/MEG file importers
% directly supported formats are CTF, Neuromag, EEP, BrainVision,
% Neuroscan and Neuralynx.
%
% Use as
%   hdr = read_fcdc_header(filename)
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
% See also READ_FCDC_DATA, READ_FCDC_EVENT

% Copyright (C) 2003-2006, Robert Oostenveld, F.C. Donders Centre
%
% $Log: read_fcdc_header.m,v $
% Revision 1.48  2009/02/17 11:08:31  roboos
% changed varargin into filename for more strict input handling: multiple input arguments are not supported because the (old) read_fcdc_xxx function interface did not support that.
%
% Revision 1.47  2009/01/28 14:45:46  roboos
% added fieldtripdefs (thanks to Verena)
%
% Revision 1.46  2008/12/16 21:21:49  roboos
% reverted to the original interface, i.e. not with variable keyval input but with a fixed list
%
% Revision 1.45  2006/06/19 10:33:57  roboos
% replaced all functional code by a call to the new low-level "fileio" function, it is now only an empty wrapper
% updated the documentation
%

fieldtripdefs

% use the low-level reading function
[hdr] = read_header(filename);

