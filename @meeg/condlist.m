function res = condlist(this)
% Method for getting a list of unique condition labels sorted according to 
% the trial order in the file
% FORMAT res = condlist(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: condlist.m 1855 2008-06-26 11:00:01Z vladimir $

res = getset(this, 'trials', 'label');

if ~iscell(res)
    res = {res};
end

[res, ind] = unique(res);

[junk, ind] = sort(ind);

res = res(ind);
