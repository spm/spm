function display(obj)
% Display method for slice overlay object
%__________________________________________________________________________

% Copyright (C) 2005-2022 Matthew Brett


X = struct(obj);
src = '[slice overlay object]';
if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    disp(src);
    disp(X)
else
    disp(' ')
    disp([inputname(1) ' =']);
    disp(' ');
    disp(src);
    disp(' ');
    disp(X)
end