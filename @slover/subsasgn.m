function result = subsasgn(this, Struct, rhs)
% Method to overload . notation in assignments.
% . assignment works directly on object fields
%__________________________________________________________________________

% Copyright (C) 2005-2022 Matthew Brett


result = builtin('subsasgn', this, Struct, rhs);
