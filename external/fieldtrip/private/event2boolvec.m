function boolvec = event2boolvec(event, varargin)

% EVENT2BOOLVEC converts between two representations of events or trials.
%
% FieldTrip uses a number of representations for events that are conceptually very similar
%   event    = structure with type, value, sample, duration and offset
%   trl      = Nx3 numerical array with begsample, endsample, offset
%   trl      = table with 3 columns for begsample, endsample, offset
%   artifact = Nx2 numerical array with begsample, endsample
%   artifact = table with 2 columns for begsample, endsample
%   boolvec  = 1xNsamples boolean vector with a thresholded TTL/trigger sequence
%   boolvec  = MxNsamples boolean matrix with a thresholded TTL/trigger sequence
%
% If trl or artifact are represented as a MATLAB table, they can have additional
% columns. These additional columns have to be named and are not restricted to
% numerical values.
%
% See also ARTIFACT2BOOLVEC, ARTIFACT2EVENT, ARTIFACT2TRL, BOOLVEC2ARTIFACT, BOOLVEC2EVENT, BOOLVEC2TRL, EVENT2ARTIFACT, EVENT2BOOLVEC, EVENT2TRL, TRL2ARTIFACT, TRL2BOOLVEC, TRL2EVENT

% Copyright (C) 2009, Ingrid Nieuwenhuis
% Copyright (C) 2020, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% get the optional input arguments or set defaults
endsample = ft_getopt(varargin, 'endsample', []); % see below for the default
type      = ft_getopt(varargin, 'type', {}); 
value     = ft_getopt(varargin, 'value', {});

if isempty(type) && isempty(value)
  % the default is to make one row per event type
  type = unique({event.type});
elseif ~isempty(type) && ~isempty(value)
  ft_error('you should specify either type or value, but not both');
end

for i=1:numel(event)
  if isempty(event(i).duration) || event(i).duration < 1
    % set the default and the minimum duration to one sample
    event(i).duration = 1;
  end
  if isempty(event(i).offset)
    % set the default to 0
    event(i).offset = 0;
  end
end

% set the default
if isempty(endsample)
  endsample = [event.sample] + [event.duration] - 1;
  endsample = max(endsample);
end

% set the default
if isempty(value)
  if all(cellfun(@ischar, {event.value}))
    value = unique({event.value});
  else
    value = {};
  end
end

boolvec = false(length(type), endsample);

for i=1:length(type)
  sel = strcmp({event.type}, type{i});
  begsample = [event(sel).sample];
  duration  = [event(sel).duration];
  if isempty(duration)
    endsample = begsample;
  elseif all(duration==0)
    endsample = begsample;
  elseif all(duration==1)
    endsample = begsample;
  else
    endsample = begsample + duration - 1;
  end
  for j=1:length(begsample)
    boolvec(i, begsample(j):endsample(j)) = true;
  end
end

for i=1:length(value)
  sel = strcmp({event.value}, value{i});
  begsample = [event(sel).sample];
  duration  = [event(sel).duration];
  if isempty(duration)
    endsample = begsample;
  elseif all(duration==0)
    endsample = begsample;
  elseif all(duration==1)
    endsample = begsample;
  else
    endsample = begsample + duration - 1;
  end
  for j=1:length(begsample)
    boolvec(i, begsample(j):endsample(j)) = true;
  end
end
