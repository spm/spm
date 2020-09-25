function res = copy(this, newname)
% Method for copying a dataset
% FORMAT res = copy(this, fname)
%
% fname can be
% - path\filename -> data copied and renamed
% - path          -> data copied only
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: copy.m 7962 2020-09-25 12:06:47Z vladimir $

[p, f] = fileparts(newname);

if ~isempty(p)
    if ~exist(p,'dir'), spm_mkdir(p); end;
else
    p = path(this);
end

if isempty(f)
    f = fname(this);
else
    f = [f '.mat'];
end

if strcmpi(fullfile(this), fullfile(p, f))
    res = this;
    return;
end

%-Copy dataset (.mat and .dat)
%--------------------------------------------------------------------------
new = clone(this, fullfile(p, f));
[r, msg] = copyfile(fnamedat(this), ...
    fnamedat(new), 'f');
if ~r
    error(msg);
    res = [];
else
    res = new;
end
