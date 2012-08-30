function str = showdetail(item)

% function str = showdetail(item)
% Display details for a cfg_repeat and all of its options.
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
str{end+1} = 'Class  : cfg_repeat';
% Display detailed help for each cfg_choice value item
citems = subsref(item, substruct('.','values'));
if numel(citems) > 1 || item.forcestruct
    str{end+1} = ['A cell array, each cell containing a struct with a ' ...
                  'single field:'];
    for k = 1:numel(citems)
        str{end+1} = sprintf('.%s', gettag(citems{k}));
    end;
else
    if isa(citems{1}, 'cfg_branch')
        str{end+1} = 'A struct array with fields:';
    else
        str{end+1} = 'A cell array with contents:';
    end
    str = [str; showdetail(citems{1})];
end