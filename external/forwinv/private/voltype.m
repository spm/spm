function [type] = voltype(vol, desired)

% VOLTYPE determines the type of volume conduction model
%
% Use as
%   [type] = voltype(vol)
% to get a string describing the type, or
%   [flag] = voltype(vol, desired)
% to get a boolean value.
%
% See also READ_VOL, COMPUTE_LEADFIELD

% Copyright (C) 2007-2008, Robert Oostenveld
%
% $Log: voltype.m,v $
% Revision 1.3  2008/03/18 12:39:24  roboos
% change in comment, nothing functional
%
% Revision 1.2  2008/03/05 15:24:44  roboos
% changed the detection of various spherical models, added infinite vacuum
%
% Revision 1.1  2007/07/25 08:31:12  roboos
% implemented new helper function
%

if isfield(vol, 'type')
  % preferably the structure specifies its own type
  type = vol.type;

elseif isfield(vol, 'r') && numel(vol.r)==1 && ~isfield(vol, 'label')
  type = 'singlesphere';

elseif isfield(vol, 'r') && numel(vol.r)>=2 && ~isfield(vol, 'label')
  type = 'concentric';

elseif isfield(vol, 'r') && isfield(vol, 'o') && isfield(vol, 'label')
  type = 'multisphere';

elseif isfield(vol, 'bnd')
  type = 'bem';

elseif isempty(vol)
  type = 'infinite';
end

if nargin>1
  type = strcmp(type, desired);
end

