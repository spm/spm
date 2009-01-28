function [dat] = read_fcdc_data(filename, header, begsample, endsample, chanindx, continuous)

% READ_FCDC_DATA is a wrapper around different EEG/MEG file importers
% directly supported formats are CTF, Neuromag, EEP, BrainVision,
% Neuroscan and Neuralynx.
%
% Use as
%   dat = read_fcdc_data(filename, header, begsample, endsample, chanindx, continuous)
% with the following inputs
%   'header'         header structure, see READ_FCDC_HEADER
%   'begsample'      first sample to read
%   'endsample'      last sample to read
%   'chanindx'       list with channel indices to read
%   'continuous'     boolean flag, whether the data should be interpreted as continuous
%
% This function treats all data as continuous, even if it is stored
% in disconnected data segments (i.e. trials or epochs). If the data
% is stored in trials and you want to read over a trial-boundary, you
% have to specify continuous as true. This is usefull for example  in
% the case od pseudo-continuous CTF data files. This function returns
% a 2-D matrix of size Nchans*Nsamples.
%
% See also READ_FCDC_HEADER, READ_FCDC_EVENT

% Copyright (C) 2003-2006, Robert Oostenveld, F.C. Donders Centre
%
% $Log: read_fcdc_data.m,v $
% Revision 1.44  2009/01/28 14:45:46  roboos
% added fieldtripdefs (thanks to Verena)
%
% Revision 1.43  2008/12/16 21:21:50  roboos
% reverted to the original interface, i.e. not with variable keyval input but with a fixed list
%
% Revision 1.42  2006/06/19 10:33:57  roboos
% replaced all functional code by a call to the new low-level "fileio" function, it is now only an empty wrapper
% updated the documentation
%

fieldtripdefs

% set the defaults for the optional input arguments
if nargin<5
  chanindx = [];
end
if nargin<6 || isempty(continuous)
  continuous = false;
end

% use the low-level reading function
[dat] = read_data(filename, 'header', header, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', ~continuous);

