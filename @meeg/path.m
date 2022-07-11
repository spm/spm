function res = path(this, newpath)
% Method for getting/setting path
% FORMAT res = path(this, newpath)
%__________________________________________________________________________

% Stefan Kiebel
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if nargin == 1
    res = this.path;
else
    this.path = newpath;
    res       = this;
end
