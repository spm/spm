function res = fname(this, newname)
% Method for getting/setting file name
% FORMAT res = fname(this, name)
%__________________________________________________________________________

% Stefan Kiebel
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if  nargin == 1
    res = this.fname;
else
    this.fname = [spm_file(newname, 'basename') '.mat'];
    res = this;
end
