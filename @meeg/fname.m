function res = fname(obj, name)
% Method for getting/setting file name
% FORMAT res = fname(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id$

switch nargin
    case 1
        res = getfname(this);        
    case 2
        res = setfname(this, name);
    otherwise
end

function res = getfnamedat(this)
res = this.fname;

function this = setfnamedat(this, name)
this.fname = name;
