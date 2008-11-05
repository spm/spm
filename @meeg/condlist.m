function res = condlist(this, newcondlist)
% Method for getting a list of unique condition labels sorted according to
% the trial order in the file
% FORMAT res = condlist(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: condlist.m 2446 2008-11-05 16:05:14Z vladimir $

res = getset(this, 'trials', 'label');
if ~iscell(res)
    res = {res};
end

[res, ind] = unique(res);

if nargin == 1

    [junk, ind] = sort(ind);

    res = res(ind);

    if numel(res)>1 && isfield(this.other, 'condlist') &&...
            iscell(this.other.condlist) && ~isempty(this.other.condlist)
        [sel1, sel2] = spm_match_str(this.other.condlist, res);
        res = res([sel2 setdiff(1:numel(res), sel2)]);
    end
else
    if iscell(newcondlist) && isempty(setdiff(newcondlist, res))
        this.other(1).condlist = newcondlist;
    else
        error('Expecting a cell array with condition labels as input.');
    end
    res = this;
end


