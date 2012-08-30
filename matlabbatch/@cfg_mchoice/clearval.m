function item = clearval(item, dflag)

% function item = clearval(item, dflag)
% Clear val field, thereby removing the currently selected configuration
% subtree. If dflag is set, then also all val fields in the item.values{:}
% cfg_item(s) are cleared. 
% This function is identical for cfg_choice and cfg_repeat items.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: clearval.m 4867 2012-08-30 13:04:51Z volkmar $

rev = '$Rev: 4867 $'; %#ok

item.cfg_item.val = {};

if dflag % clear defaults
    for k = 1:numel(item.values)
        item.values{k} = clearval(item.values{k}, ...
            dflag);
    end;
end;
