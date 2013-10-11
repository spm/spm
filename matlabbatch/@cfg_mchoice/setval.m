function item = setval(item, val, dflag)

% function item = setval(item, val, dflag)
% If isempty(val), set item.val to {}. Otherwise, if item.values{val(1)}
% is not already in item.val, set item.val{end+1} to item.values{val(1)}.
% This method does not yet delete individual val items.
% dflag is ignored for cfg_mchoice items.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: setval.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

if isempty(val)
    item = subsasgn(item, substruct('.','val'), {});
else
    if ~any(strcmp(gettag(item.values{val(1)}), tagnames(item, false)))
        val = item.values{val(1)};
        item = subsasgn(item, substruct('.','val', '{}',{numel(item.cfg_item.val)+1}), val);
    end;
end;
