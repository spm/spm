function ok = all_set_item(item)

% function ok = all_set_item(item)
% Perform within-item all_set check. For repeats, this is true, if item.val
% has between item.num(1) and item.num(2) elements.
%
% This code is part of a batch job configuration system for MATLAB. See
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: all_set_item.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok
ok = (numel(item.cfg_item.val) >= item.num(1)) &&...
    (numel(item.cfg_item.val) <= item.num(2));
