function [r] = issubfield(s, f);

% ISSUBFIELD tests for the presence of a field in a structure just like the standard
% Matlab ISFIELD function, except that you can also specify nested fields
% using a '.' in the fieldname. The nesting can be arbitrary deep.
%
% Use as
%   f = issubfield(s, 'fieldname')
% or as
%   f = issubfield(s, 'fieldname.subfieldname')
%
% This function returns true if the field is present and false if the field
% is not present.
%
% See also ISFIELD, GETSUBFIELD, SETSUBFIELD

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: issubfield.m,v $
% Revision 1.1  2005/02/08 09:30:18  roboos
% new implementations to make it easier to work with nested structures
%

try
  getsubfield(s, f);    % if this works, then the subfield must be present  
  r = 1;
catch
  r = 0;                % apparently the subfield is not present
end