function [sens] = apply_montage(sens, montage, varargin)

% APPLY_MONTAGE changes the montage of an electrode or gradiometer array. A
% montage can be used for EEG rereferencing, MEG synthetic gradients, MEG
% planar gradients or unmixing using ICA. This function applies the montage
% to the sensor array. The sensor array can subsequently be used for
% forward computation and source reconstruction of the data.
%
% Use as
%   [sens] = apply_montage(sens, montage, ...)
%
% A montage is specified as a structure with the fields
%   montage.tra      = MxN matrix
%   montage.labelnew = Mx1 cell-array
%   montage.labelorg = Nx1 cell-array
%
% Additional options should be specified in key-value pairs and can be
%   'keepunused'   string, 'yes' or 'no' (default = 'no')
%   'inverse'      string, 'yes' or 'no' (default = 'no')
%
% See also READ_SENS, TRANSFORM_SENS

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: apply_montage.m,v $
% Revision 1.5  2008/05/15 15:08:51  roboos
% added support for applying the inverse montage (i.e. undo a previous montage)
% added support for applying the montage to preprocessed/raw data
%
% Revision 1.4  2008/05/13 11:43:27  roboos
% fixed bug in selempty
%
% Revision 1.3  2008/05/13 09:08:19  roboos
% fixed bug in assignment
% added option keepunused=yes|no
%
% Revision 1.2  2008/04/28 14:14:29  roboos
% added closing bracket
%
% Revision 1.1  2008/04/14 20:16:37  roboos
% new implementation, required for spm integration
%

% get optional input arguments
keepunused = keyval('keepunused', varargin{:}); if isempty(keepunused), keepunused = 'no'; end
inverse    = keyval('inverse',    varargin{:}); if isempty(inverse),    inverse    = 'no'; end

if strcmp(inverse, 'yes')
  % apply the inverse montage, i.e. undo a previously applied montage
  tmp.labelnew = montage.labelorg;
  tmp.labelorg = montage.labelnew;
  tmp.tra      = inv(montage.tra);
  montage      = tmp;
end

% use default transfer from sensors to channels if not specified
if ~isfield(sens, 'tra') && ~isfield(sens, 'trial')
  nchan = size(sens.pnt,1);
  sens.tra = sparse(eye(nchan));
end

% select and discard the columns that are empty
selnonempty       = find(~all(montage.tra==0, 1));
montage.tra       = montage.tra(:,selnonempty);
montage.labelorg  = montage.labelorg(selnonempty);

% add columns for the channels that are not involved in the montage
add = setdiff(sens.label, montage.labelorg);
m = size(montage.tra,1);
n = size(montage.tra,2);
k = length(add);
if strcmp(keepunused, 'yes')
  % add the channels that are not rereferenced to the input and output
  montage.tra((m+(1:k)),(n+(1:k))) = eye(k);
  montage.labelorg = cat(1, montage.labelorg(:), add(:));
  montage.labelnew = cat(1, montage.labelnew(:), add(:));
else
  % add the channels that are not rereferenced to the input montage only
  montage.tra(:,(n+(1:k))) = zeros(m,k);
  montage.labelorg = cat(1, montage.labelorg(:), add(:));
end

% determine whether all channels are unique
m = size(montage.tra,1);
n = size(montage.tra,2);
if length(unique(montage.labelnew))~=m
  error('not all output channels of the montage are unique');
end
if length(unique(montage.labelorg))~=n
  error('not all input channels of the montage are unique');
end

% determine whether all channels that have to be rereferenced are available
if length(intersect(sens.label, montage.labelorg))~=length(sens.label)
  error('not all channels that are used in the montage are available');
end

% reorder the columns of the montage matrix
[selsens, selmont] = match_str(sens.label, montage.labelorg);
montage.tra        = sparse(montage.tra(:,selmont));
montage.labelorg   = montage.labelorg(selmont);

if isfield(sens, 'tra')
  % apply the montage to the sensor array
  sens.tra   = montage.tra * sens.tra;
  sens.label = montage.labelnew;
elseif isfield(sens, 'trial')
  % apply the montage to the data that was preprocessed using fieldtrip
  Ntrials = numel(sens.trial);
  for i=1:Ntrials
    fprintf('processing trial %d from %d\n', i, Ntrials);
    sens.trial{i}   = montage.tra * sens.trial{i};
  end
  sens.label = montage.labelnew;
else
  error('unrecognized input');
end
