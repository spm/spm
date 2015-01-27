function res = path(this, newpath)
% Method for getting/setting path
% FORMAT res = path(this, newpath)
% _______________________________________________________________________
% Copyright (C) 2008-2015 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: path.m 6318 2015-01-27 10:36:34Z vladimir $

if nargin == 1
    res = this.path;
else
    this.path = newpath;
    res       = this;
end

