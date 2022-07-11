function result = subsref(this, Struct)
% Method to overload the . notation.
% . reference works directly on object fields
%__________________________________________________________________________

% Copyright (C) 2005-2022 Matthew Brett


result = builtin('subsref', this, Struct );
