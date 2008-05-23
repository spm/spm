function str = showdoc(item, indent)

% function str = showdoc(item, indent)
% Display help text for a cfg_repeat and all of its options.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: showdoc.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

str = showdoc(item.cfg_item, indent);
citems = subsref(item, substruct('.','values'));
if isfinite(item.num(2))
    str{end+1} = sprintf('Between %d and %d options must be selected from:', ...
                         item.num(1), item.num(2));
else
    str{end+1} = sprintf('%d or more options must be selected from:', ...
                         item.num(1));
end;
str{end+1} = '';
% Display short listing of cfg_repeat value items first
for k = 1:numel(citems)
    str{end+1} = sprintf('* %s', subsref(citems{k}, substruct('.','name')));
end;
valitem = subsref(item, substruct('.','val'));
if ~isempty(valitem)
    str{end+1} = sprintf('Currently selected options:');
    for k = 1:numel(valitem)
        str{end+1} = sprintf('* "%s"', subsref(valitem{k}, substruct('.','name')));
    end;
end;
str{end+1} = '';
% Display detailed help for each cfg_choice value item
for k = 1:numel(citems)
    str1 = showdoc(citems{k}, sprintf('%s%d.', indent, k));
    str = {str{:} str1{:}};
    str{end+1} = '';
end;