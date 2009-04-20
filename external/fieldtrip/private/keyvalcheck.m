function keyvalcheck(arglist, varargin)

% KEYVALCHECK is a helper function for parsing optional key-value input pairs.
%
% Use as
%   keyvalcheck(argin, 'required',  {'key1', 'key2', ...})
%   keyvalcheck(argin, 'forbidden', {'key1', 'key2', ...})
%   keyvalcheck(argin, 'optional',  {'key1', 'key2', ...})
%
% See also KEYVAL

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: keyvalcheck.m,v $
% Revision 1.1  2009/04/14 19:37:44  roboos
% new helper function, used in plot_xxx functions
%

required  = keyval('required', varargin);
forbidden = keyval('forbidden', varargin);
optional  = keyval('optional', varargin);

keys = arglist(1:2:end);
vals = arglist(2:2:end);

if numel(keys)~=numel(vals)
  error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
end

set = intersect(keys, required);
if numel(set)~=numel(required)
  error('the required input argument ''%s'' was not specified', set{:});
end

set = intersect(keys, forbidden);
if numel(set)~=0
  error('the input argument ''%s'' is forbidden', set{:});
end

set = setdiff(keys, optional);
if numel(set)>0
  error('the input argument ''%s'' is forbidden', set{:});
end
