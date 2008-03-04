function item = setval(item, val)

% function item = setval(item, val)
% Set item.val{1} to item.values{val}. If val == {}, set item.val to {}.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: setval.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

if iscell(val) && isempty(val)
    item.val = {};
else
    val = item.values{val};
    item = subsasgn(item, substruct('.','val', '{}',{1}), val);
end;
