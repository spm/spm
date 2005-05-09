function result = subsref(this, Struct)
% method to overload the . notation.
% . reference works directly on object fields
%
% $Id: subsref.m,v 1.1 2005/04/20 15:05:36 matthewbrett Exp $

result = builtin('subsref', this, Struct );