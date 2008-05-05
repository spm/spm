function str = showdoc(item, indent)

% function str = showdoc(item, indent)
% Display help text for a cfg_branch and all of its children.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: showdoc.m 1541 2008-05-05 13:36:51Z volkmar $

rev = '$Rev: 1541 $';

str = showdoc(item.cfg_item, indent);
citems = subsref(item, substruct('.','val'));
str{end+1} = '';
str{end+1} = sprintf('This branch contains %d items:', numel(citems));
% Display short listing of branch items first
for k = 1:numel(citems)
    str{end+1} = sprintf('* %s', subsref(citems{k}, substruct('.','name')));
end;
str{end+1} = '';
% Display detailed help for each branch item
for k = 1:numel(citems)
    str1 = showdoc(citems{k}, sprintf('%s%d.', indent, k));
    str = {str{:} str1{:}};
    str{end+1} = '';
end;