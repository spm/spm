function res = fnamedat(this)
% Method for getting the name of the data file
% FORMAT res = fnamedat(this)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: fnamedat.m 5025 2012-10-31 14:44:13Z vladimir $


if islinked(this)
    res = this.data.fname;
else
    res = [];
end