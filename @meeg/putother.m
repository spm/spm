function obj = putother(obj, stuff)
% Method for putting a variable into the other field
% FORMAT res = putother(obj, stuff)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

eval(['obj.other(1).' inputname(2) ' = stuff;']);
