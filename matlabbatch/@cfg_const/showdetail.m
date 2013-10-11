function str = showdetail(item)

% function str = showdetail(item)
% Display details for a cfg_const item.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: showdetail.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

str = showdetail(item.cfg_item);
str{end+1} = 'This item has a constant value which can not be modified.';
