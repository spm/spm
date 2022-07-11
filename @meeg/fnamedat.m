function res = fnamedat(this)
% Method for getting the name of the data file
% FORMAT res = fnamedat(this)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if islinked(this)
    res = this.data.fname;
else
    res = [];
end