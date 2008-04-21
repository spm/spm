function tree = setfilename(tree,filename)
% XMLTREE/SETFILENAME Set filename method
% FORMAT tree = setfilename(tree,filename)
% 
% tree     - XMLTree object
% filename - XML filename
%_______________________________________________________________________
%
% Set the filename linked to the XML tree as filename.
%_______________________________________________________________________
% Copyright (C) 2002-2008  http://www.artefact.tk/

% Guillaume Flandin <guillaume@artefact.tk>
% $Id: setfilename.m 1460 2008-04-21 17:43:18Z guillaume $

tree.filename = filename;
