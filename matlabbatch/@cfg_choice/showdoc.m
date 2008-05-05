function str = showdoc(item, indent)

% function str = showdoc(item, indent)
% Display help text for a cfg_choice and all of its options.
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
str{end+1} = '';
citems = subsref(item, substruct('.','values'));
str{end+1} = sprintf('One of the following options must be selected:');
% Display short listing of cfg_choice value items first
for k = 1:numel(citems)
    str{end+1} = sprintf('* %s', subsref(citems{k}, substruct('.','name')));
end;
valitem = subsref(item, substruct('.','val'));
if ~isempty(valitem)
    str{end+1} = sprintf('Currently selected option:');
    str{end+1} = sprintf('* "%s"',subsref(valitem{1}, substruct('.','name')));
end;
str{end+1} = '';
% Display detailed help for each cfg_choice value item
for k = 1:numel(citems)
    str1 = showdoc(citems{k}, sprintf('%s%d.', indent, k));
    str = {str{:} str1{:}};
    str{end+1} = '';
end;