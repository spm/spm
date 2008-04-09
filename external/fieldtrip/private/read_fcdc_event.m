function [event] = read_fcdc_event(varargin)

% READ_FCDC_EVENT reads all events from an EEG/MEG dataset and returns
% them in a well defined structure. It is a wrapper around different
% EEG/MEG file importers, directly supported formats are CTF, Neuromag,
% EEP, BrainVision, Neuroscan and Neuralynx.
%
% Use as
%   [event] = read_fcdc_event(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'eventformat'   string
%
% This function returns a structure with the following fields
%   event.type     = string
%   event.sample   = expressed in samples, first sample of file is 1
%   event.value    = number or string
%   event.offset   = expressed in samples
%   event.duration = expressed in samples
%
% Type and sample fields are always defined, other fields can be empty,
% depending on the type of event file. Events are sorted by the sample on
% which they occur. After reading the event structure, you can use the
% following tricks to extract information about those events in which you
% are interested.
%
% Determine the different event types
%   unique({event.type})
%
% Get the index of all trial events
%   find(strcmp('trial', {event.type}))
%
% Make a vector with all triggers that occurred on the backpanel
%   [event(find(strcmp('backpanel trigger', {event.type}))).value]
%
% Find the events that occurred in trial 26
%   t=26; samples_trials = [event(find(strcmp('trial', {event.type}))).sample];
%   find([event.sample]>samples_trials(t) & [event.sample]<samples_trials(t+1))
%
% See also READ_FCDC_HEADER, READ_FCDC_DATA

% Copyright (C) 2004-2006, Robert Oostenveld, F.C. Donders Centre
%
% $Log: read_fcdc_event.m,v $
% Revision 1.46  2006/06/19 10:33:57  roboos
% replaced all functional code by a call to the new low-level "fileio" function, it is now only an empty wrapper
% updated the documentation
%

% use the low-level reading function
[event] = read_event(varargin{:});

