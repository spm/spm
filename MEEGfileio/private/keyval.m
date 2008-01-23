function [val] = keyval(key, varargin);

% KEYVAL returns the value that corresponds to the requested key in a
% key-value pair list of variable input arguments
%
% Use as
%   [val] = keyval(key, varargin)
%
% See also VARARGIN

% Copyright (C) 2005-2007, Robert Oostenveld
%
% $Log: keyval.m,v $
% Revision 1.2  2007/07/18 12:43:53  roboos
% test for an even number of optional input arguments
%
% Revision 1.1  2005/11/04 10:24:46  roboos
% new implementation
%

if length(varargin)==1 && iscell(varargin{1})
  varargin = varargin{1};
end

if mod(length(varargin),2)
  error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
end

keys = varargin(1:2:end);
vals = varargin(2:2:end);

hit = find(strcmp(key, keys));
if length(hit)==0
  % the requested key was not found
  val = [];
elseif length(hit)==1  
  % the requested key was  found
  val = vals{hit};
else
  error('multiple input arguments with the same name');
end

