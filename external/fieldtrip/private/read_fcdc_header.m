function [hdr] = read_fcdc_header(varargin);

% READ_FCDC_HEADER is a wrapper around different EEG/MEG file importers
% directly supported formats are CTF, Neuromag, EEP, BrainVision,
% Neuroscan and Neuralynx.
%
% Use as
%   hdr = read_fcdc_header(filename, ...)
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
% See also READ_FCDC_DATA, READ_FCDC_EVENT

% Copyright (C) 2003-2006, Robert Oostenveld, F.C. Donders Centre
%
% $Log: read_fcdc_header.m,v $
% Revision 1.45  2006/06/19 10:33:57  roboos
% replaced all functional code by a call to the new low-level "fileio" function, it is now only an empty wrapper
% updated the documentation
%

% use the low-level reading function
[hdr] = read_header(varargin{:});

