function res = path(obj, name)
% Method for getting/setting path
% FORMAT res = path(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id$

switch nargin
    case 1
        res = getpath(this);        
    case 2
        res = setpath(this, name);
    otherwise
end

function res = getpath(this)
try
    res = obj.path;
catch
    res = '';
end

function this = setpath(this, name)
this.path = name;


