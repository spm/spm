function res = transformtype(this, newtype)
% Method for getting/setting type of transform
% FORMAT res = transformtype(this, name)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: transformtype.m 2081 2008-09-11 13:04:24Z vladimir $

if nargin == 1
    res = this.transform.ID;
else
    this.transform.ID = newtype;
    res = this;
end
