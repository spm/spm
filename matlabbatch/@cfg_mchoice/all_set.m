function ok = all_set(item)

% function ok = all_set(item)
% Return true, if all child items in item.val{:} are set and item specific
% criteria (i.e. number of element in .val) are met. No checks based on
% the content of item.val are performed here. 
% Content checking is done in the following places:
% * context-insensitive checks based on configuration specifications
%   are performed during subsasgn/setval. This will happen during user
%   input or while resolving dependencies during harvest. 
% * context sensitive checks by a configuration .check function are
%   performed during harvest after all dependencies are resolved.
% This function is identical for all in-tree items.
%
% This code is part of a batch job configuration system for MATLAB. See
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: all_set.m 4867 2012-08-30 13:04:51Z volkmar $

rev = '$Rev: 4867 $'; %#ok

ok = all_set_item(item);
if ok
    for k = 1:numel(item.cfg_item.val)
        ok = ok && all_set(item.cfg_item.val{k});
        if ~ok
            break;
        end;
    end;
end;