function str = showdoc(item, indent)

% function str = showdoc(item, indent)
% Generic showdoc function for cfg_item classes. It displays the
% (indented) name of the item and the justified help text for this item.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: showdoc.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

str = {sprintf('%s %s', indent, item.name)};
if ~isempty(item.help)
    str = {str{:} item.help{:}};
end;
if ~isempty(item.def)
    str{end+1} = sprintf(['This item has a default value, set via a call ' ...
                        'to function']);
    argstr = sprintf('arg%d,', 1:numel(item.def)-1);
    str{end+1} = sprintf('   %s(%s)',func2str(item.def{1}),argstr(1:end-1));
    str{end+1} = 'with arguments';
    for k = 2:numel(item.def)
        str1 = gencode(item.def{k},sprintf('   arg%d',k-1));
        str = {str{:} str1{:}};
    end;
end;