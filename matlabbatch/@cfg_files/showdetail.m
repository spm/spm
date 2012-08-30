function str = showdetail(item)

% function str = showdetail(item)
% Display details for a cfg_files item.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: showdetail.m 4867 2012-08-30 13:04:51Z volkmar $

rev = '$Rev: 4867 $'; %#ok

str = showdetail(item.cfg_item);
str{end+1} = 'Class  : cfg_files';
str{end+1} = 'A cellstr array of file names.';
str = [str; gencode(item.filter, 'filter:')];
str = [str; gencode(item.num,    'num   :')];