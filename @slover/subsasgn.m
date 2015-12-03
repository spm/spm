function result = subsasgn(this, Struct, rhs)
% Method to overload . notation in assignments.
% . assignment works directly on object fields
%__________________________________________________________________________

% Matthew Brett
% $Id: subsasgn.m 6623 2015-12-03 18:38:08Z guillaume $

result = builtin('subsasgn', this, Struct, rhs);
