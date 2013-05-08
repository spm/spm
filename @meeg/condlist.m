function res = condlist(this, newcondlist)
% Method for getting a list of unique condition labels sorted according to
% the trial order in the file
% FORMAT res = condlist(this)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: condlist.m 5472 2013-05-08 00:24:36Z vladimir $

res = getset(this, 'trials', 'label');

if isempty(res)
    if nargin == 1
        res = {};
    else
        res = this;
    end
    return;
end
    

if ~iscell(res)
    res = {res};
end


[res, ind] = unique(res);


if nargin == 1

    [junk, ind] = sort(ind);

    res = res(ind);

    if numel(res)>1  && ~isempty(this.condlist)
        [sel1, sel2] = match_str(this.condlist, res);
        res = res([sel2(:)' setdiff(1:numel(res), sel2)]);
    end
else
    if iscell(newcondlist) && ~isempty(intersect(newcondlist, res))  
        [junk, ind] = unique(newcondlist, 'first');
        newcondlist = newcondlist(sort(ind));
        
        this.condlist = newcondlist(ismember(newcondlist, res));
    else
        error('Expecting a cell array with condition labels as input.');
    end
    res = this;
end


