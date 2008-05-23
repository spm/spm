function ok = all_leafs(item)

% function ok = all_leafs(item)
% Return true, if all child items in item.val{:} consist of subtrees
% ending in leaf nodes. Leaf nodes do not have to be set at this time and
% no checks on their contents will be performed.
% This function is identical for all in-tree items. 
%
% This code is part of a batch job configuration system for MATLAB. See
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: all_leafs.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

% Check whether the required number of val fields is present. This is
% done in all_set_item for in-tree nodes.
ok = all_set_item(item);
if ok
    for k = 1:numel(item.cfg_item.val)
        ok = ok && all_leafs(item.cfg_item.val{k});
        if ~ok
            break;
        end;
    end;
end;
