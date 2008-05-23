function item = setval(item, val, dflag)

% function item = setval(item, val, dflag)
% Set item.val{1} to item.values{val(1)}. If isempty(val), set item.val to {}.
% dflag is ignored for cfg_choice items.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: setval.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

if isempty(val)
    item = subsasgn(item, substruct('.','val'), {});
else
    val = item.values{val(1)};
    item = subsasgn(item, substruct('.','val', '{}',{1}), val);
end;
