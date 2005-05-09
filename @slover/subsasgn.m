function result = subsasgn(this, Struct, rhs)
% method to overload . notation in assignments.
% . assignment works directly on object fields
%
% $Id: subsasgn.m,v 1.1 2005/04/20 15:05:36 matthewbrett Exp $

result = builtin('subsasgn', this, Struct, rhs);
