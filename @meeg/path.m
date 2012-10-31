function res = path(this, newpath)
% Method for getting/setting path
% FORMAT res = path(this, newpath)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: path.m 5025 2012-10-31 14:44:13Z vladimir $

if nargin == 1
    res = this.path;
else
    this.path = newpath;
end

