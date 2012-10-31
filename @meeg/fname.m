function res = fname(this, newname)
% Method for getting/setting file name
% FORMAT res = fname(this, name)
% _______________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: fname.m 5025 2012-10-31 14:44:13Z vladimir $

if  nargin == 1
    res = this.fname;
else
    this.fname = [spm_file(newname, 'basename') '.mat'];
    res = this;
end
