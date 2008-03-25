function res = fnamedat(obj, name)
% Method for getting/setting file name of data file
% FORMAT res = fnamedat(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id$


switch nargin
    case 1
        res = getfnamedat(this);        
    case 2
        res = setfnamedat(this, name);
    otherwise
end

function res = getfnamedat(this)
res = obj.data.fnamedat;

function this = setfnamedat(this, name)
this.data.fnamedat = fnamedat;
