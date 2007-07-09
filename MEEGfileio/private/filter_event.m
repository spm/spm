function event = filter_event(event, varargin)

% FILTER_EVENT does what its name implies
%
% Use as
%   event = filter_event(event, ...)
%
% The optional arguments should come in key-value pairs and determine the
% filter characteristics.
%
% See also READ_EVEN, WRITE_EVENT

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: filter_event.m,v $
% Revision 1.1  2007/06/13 14:47:35  roboos
% moved filter_event from fieldtrip to fileio module
%
% Revision 1.1  2007/06/06 12:41:22  roboos
% new implementation
%

% get the optional input arguments
type      = keyval('type', varargin);
value     = keyval('value', varargin);
sample    = keyval('sample', varargin);
timestamp = keyval('timestamp', varargin);
offset    = keyval('offset', varargin);
duration  = keyval('duration', varargin);

% the numeric fields can also be filtered on a range
minsample = keyval('minsample', varargin);
maxsample = keyval('maxsample', varargin);
minduration = keyval('minduration', varargin);
maxduration = keyval('maxduration', varargin);

if ~isempty(type)
  % this can be specified as string or as cell-array, convert to cell-array
  if ~iscell(type)
    type = {type};
  end
end

% determine which filters to apply
testtype        = ~isempty(type) && isfield(event, 'type');
testvalue       = ~isempty(value) && isfield(event, 'value');
testsample      = ~isempty(sample) && isfield(event, 'sample');
testtimestamp   = ~isempty(timestamp) && isfield(event, 'timestamp');
testoffset      = ~isempty(offset) && isfield(event, 'offset');
testduration    = ~isempty(duration) && isfield(event, 'duration');
testminsample   = ~isempty(minsample) && isfield(event, 'sample');
testmaxsample   = ~isempty(maxsample) && isfield(event, 'sample');
testminduration = ~isempty(minduration) && isfield(event, 'duration');
testmaxduration = ~isempty(maxduration) && isfield(event, 'duration');
% TODO mintimestamp, maxtimestamp

% apply the filters
sel = true(length(event),1);
for i=1:length(event)
  % test whether they match with the selected arrays
  if testsample,   sel(i) = sel(i) && any(event(i).sample == sample); end
  if testoffset,   sel(i) = sel(i) && any(event(i).offset == offset); end
  if testduration, sel(i) = sel(i) && any(event(i).duration == duration); end
  % test whether they ly within the specified range
  if testminsample,   sel(i) = sel(i) && any(event(i).sample>=minsample); end
  if testmaxsample,   sel(i) = sel(i) && any(event(i).sample>=maxsample); end
  if testminduration, sel(i) = sel(i) && any(event(i).sample>=minduration); end
  if testmaxduration, sel(i) = sel(i) && any(event(i).sample>=maxduration); end
  % this is potentially the slowest test, hence do it the last
  if testtype,     sel(i) = sel(i) && any(strcmp(event(i).type, type)); end
end

event = event(sel);
