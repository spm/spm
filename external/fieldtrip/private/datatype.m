function [type, dimord] = datatype(data, desired)

% DATATYPE determines the type of data represented in a FieldTrip data
% structure and returns a string with raw, freq, timelock source, comp,
% spike, source, volume, dip.
%
% Use as
%   [type, dimord] = datatype(data)
%   [type, dimord] = datatype(data, desired)
%
% See also CHANTYPE, FILETYPE, SENSTYPE, VOLTYPE

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: datatype.m,v $
% Revision 1.2  2008/12/19 09:12:54  roboos
% added support for desired type, returning boolean
%
% Revision 1.1  2008/12/18 15:49:55  roboos
% new implementation
%

% determine the type of input data, this can be raw, freq, timelock, comp, spike, source, volume, dip
israw      = isfield(data, 'label') && isfield(data, 'time') && isa(data.time, 'cell') && isfield(data, 'trial') && isa(data.trial, 'cell');
isfreq     = isfield(data, 'label') && isfield(data, 'freq') && (isfield(data, 'powspctrm') || isfield(data, 'crsspctrm') || isfield(data, 'fourierspctrm'));
istimelock = isfield(data, 'label') && isfield(data, 'time') && ~isfield(data, 'freq') && isfield(data, 'avg') && isnumeric(data.avg);
iscomp     = isfield(data, 'topo') || isfield(data, 'topolabel');
isspike    = isfield(data, 'label') && isfield(data, 'waveform') && isa(data.waveform, 'cell') && isfield(data, 'timestamp') && isa(data.timestamp, 'cell');
isvolume   = isfield(data, 'transform') && isfield(data, 'dim');
issource   = isfield(data, 'pos');
isdip      = isfield(data, 'dip');

if israw
  type = 'raw';
elseif isfreq
  type = 'freq';
elseif istimelock
  type = 'timelock';
elseif iscomp
  type = 'comp';
elseif isspike
  type = 'spike';
elseif isvolume
  type = 'volume';
elseif issource
  type = 'source';
elseif isdip
  type = 'dip';
end

if nargin>1
  % return a boolean value
  type = strcmp(type, desired);
  return;
end

if nargout>1
  % also return the dimord of the input data
  if isfield(data, 'dimord')
    dimord = data.dimord;
  else
    dimord = 'unknown';
  end
end

