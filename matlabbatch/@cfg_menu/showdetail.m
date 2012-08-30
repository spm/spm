function str = showdetail(item)

% function str = showdetail(item)
% Display details for a cfg_menu and all of its options.
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
str{end+1} = 'Class  : cfg_menu';
str{end+1} = sprintf('The value for one of the following options must be entered:');
% Display cfg_menu labels
for k = 1:numel(item.labels)
    str(end+1) = gencode(item.values{k}, sprintf('''%s:''', item.labels{k}));
end;
