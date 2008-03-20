function obj = putcache(obj, stuff)
% Method for putting a variable into a temporary cache
% FORMAT res = putcache(obj, stuff)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

eval(['obj.cache(1).' inputname(2) ' = stuff;']);
