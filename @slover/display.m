function display(obj)
% Display method for slice overlay object
%__________________________________________________________________________

% Matthew Brett
% $Id: display.m 6623 2015-12-03 18:38:08Z guillaume $

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