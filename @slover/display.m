function display(obj)
% display method for slice overlay object
%
% $Id: display.m,v 1.1 2005/04/20 15:05:36 matthewbrett Exp $
  
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