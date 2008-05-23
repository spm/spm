function [id, stop, rtaglist] = tag2cfgsubs(item, taglist, finalspec, tropts)

% function [id, stop, rtaglist] = tag2cfgsubs(item, taglist, finalspec, tropts)
% Return the index into the values branch of a configuration tree which
% corresponds to a list of tags. 
% This is the tag2cfgsubs function for within-tree cfg_items (cfg_[ex]branch,
% cfg_repeat, cfg_choice).
% Traversal stops if taglist contains only one element or item matches a
% non-empty tropts.stopspec. In this case, stop returns the match status.
% Id is an empty substruct, if gettag(item) matches taglist{1} and item
% matches finalspec, otherwise it is an empty cell.
% If taglist contains more than one element and taglist{2} matches any tag
% of a .val element, then the subscript index to this element is returned.
% If the recursive match was unsuccessful, it returns an empty cell and
% stop = true.
% rtaglist contains the remaining tags that were not matched due to a
% stopping criterion.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: tag2cfgsubs.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok
if numel(taglist) == 1 || (~isempty(tropts.stopspec) ...
                           && match(item, tropts.stopspec))
    if strcmp(gettag(item), taglist{1}) && match(item, finalspec)
        id = struct('type', {}, 'subs', {});
        stop = ~isempty(tropts.stopspec) ...
               && match(item, tropts.stopspec);
    else
        % tag does not match but traversal should stop
        id = {};
        stop = true;
    end;
    rtaglist = taglist(2:end);
    return;
end;

% search defaults tree
tname = treepart(item, true);
citems = subsref(item, substruct('.',tname));
id = {};
stop = true;
rtaglist = taglist(2:end);
for k = 1:numel(citems)
    if strcmp(gettag(citems{k}), taglist{2})
        [id stop rtaglist] = tag2cfgsubs(citems{k}, taglist(2:end), ...
                                         finalspec, tropts);
        if isstruct(id)
            id = [substruct('.', tname, '{}', {k}) id];
        end;
        % stop after first matching tag
        break;
    end;
end;