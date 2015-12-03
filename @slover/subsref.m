function result = subsref(this, Struct)
% Method to overload the . notation.
% . reference works directly on object fields
%__________________________________________________________________________

% Matthew Brett
% $Id: subsref.m 6623 2015-12-03 18:38:08Z guillaume $

result = builtin('subsref', this, Struct );
