function ok = all_leafs(item)

% function ok = all_leafs(item)
% Generic all_leafs function that returns true. This is suitable for all
% leaf items. No content specific checks are performed.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: all_leafs.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';
% do not check anything else than item class
ok = isa(item, 'cfg_leaf');

