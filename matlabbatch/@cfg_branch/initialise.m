function item = initialise(item, val, dflag)

% function item = initialise(item, val, dflag)
% Disassemble val struct and pass its fields on to a item.val{} item whose
% tag matches the current field name. dflag is ignored in a cfg_branch.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: initialise.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

% Determine possible tags
vtags = fieldnames(val);

for k = 1:numel(item.cfg_item.val)
    % find field in val that corresponds to one of the branch vals
    vi = strcmp(gettag(item.cfg_item.val{k}), vtags);
    if any(vi) % field names are unique, so there will be at most one match
        item.cfg_item.val{k} = initialise(item.cfg_item.val{k}, ...
            val.(vtags{vi}), dflag);
    end;
end;

