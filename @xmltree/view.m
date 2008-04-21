function view(tree)
% XMLTREE/VIEW View Method
% FORMAT view(tree)
% 
% tree   - XMLTree object
%_______________________________________________________________________
%
% Display an XML tree in a graphical interface
%_______________________________________________________________________
% Copyright (C) 2002-2008  http://www.artefact.tk/

% Guillaume Flandin <guillaume@artefact.tk>
% $Id: view.m 1460 2008-04-21 17:43:18Z guillaume $

error(nargchk(1,1,nargin));

editor(tree);
