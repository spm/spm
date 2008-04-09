function [dat] = read_fcdc_data(varargin);

% READ_FCDC_DATA is a wrapper around different EEG/MEG file importers
% directly supported formats are CTF, Neuromag, EEP, BrainVision,
% Neuroscan and Neuralynx.
%
% Use as
%   dat = read_data(filename, ...) or
%   dat = read_data(filename, header, begsample, endsample, chanindx, continuous)
%
% If you use the first function call, the options should be specified in key-value pairs and can be
%   'header'         header structure, see READ_FCDC_HEADER
%   'begsample'      first sample to read
%   'endsample'      last sample to read
%   'begtrial'       first trial to read, mutually exclusive with begsample+endsample
%   'endtrial'       last trial to read, mutually exclusive with begsample+endsample
%   'chanindx'       list with channel indices to read
%   'checkboundary'  boolean, whether to check for reading segments over a trial boundary
%   'dataformat'     string
%   'headerformat'   string
%
% If you use the second function call, the function treats all data
% as continuous, even if it is stored in disconnected data segments
% (i.e. trials or epochs). If you want to read the second trial, you
% should read from sample hdr.nSamples+1 up to 2*hdr.nSamples. The
% argument "chanindx" is optional, the default is to read all channels.
% In that case the argument "continuous" is optional, the default is
% check to whether the data is continuous and, if not, check whether
% the requested data segment does not extend over a trial boundary.
% If continuous=1, this check is not performed (usefull for
% pseudo-continuous data).
%
% This function returns a 2-D matrix of size Nchans*Nsamples for
% continuous data when begevent and endevent are specified, or a 3-D
% matrix of size Nchans*Nsamples*Ntrials for epoched or trial-based
% data when begtrial and endtrial are specified.
% 
% See also READ_FCDC_HEADER, READ_FCDC_EVENT

% Copyright (C) 2003-2006, Robert Oostenveld, F.C. Donders Centre
%
% $Log: read_fcdc_data.m,v $
% Revision 1.42  2006/06/19 10:33:57  roboos
% replaced all functional code by a call to the new low-level "fileio" function, it is now only an empty wrapper
% updated the documentation
%

% use the low-level reading function
[dat] = read_data(varargin{:});

