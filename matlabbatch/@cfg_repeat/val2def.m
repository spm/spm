function [item, defaults] = val2def(item, defaults, funname, deftag)
% function [item, defaults] = val2def(item, defaults, funname, deftag)
% If a cfg_leaf item has a value, extract it and generate code for defaults
% retrieval. This function works in a way similar to harvest, but with a
% much simpler logic. Also, it modifies the returned configuration tree by
% clearing the .val fields if they are moved to defaults.
% Initially, defaults and deftag should be empty.
% This function is identical for cfg_branch and cfg_choice classes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: val2def.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

csubs = substruct('.', treepart(item, true));
citems = subsref(item, csubs);
if numel(citems) > 1
    if ~isempty(deftag)
        deftag = [deftag '.'];
    end
    for k = 1:numel(citems)
        [citems{k}, defaults] = val2def(citems{k}, defaults, funname, [deftag gettag(citems{k})]);
    end
else
    [citems{1}, defaults] = val2def(citems{1}, defaults, funname, deftag);
end
item = subsasgn(item, csubs, citems);
