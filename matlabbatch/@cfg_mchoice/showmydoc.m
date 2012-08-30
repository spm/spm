function str = showmydoc(item, indent)

% function str = showmydoc(item, indent)
% Display help text for a cfg_choice and all of its options, without
% recursive calls to child nodes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: showmydoc.m 4867 2012-08-30 13:04:51Z volkmar $

rev = '$Rev: 4867 $'; %#ok
str1 = showdoc(item.cfg_item, indent);
citems = subsref(item, substruct('.','values'));
str2{1} = 'Some of the following options can be selected:';
% Display short listing of cfg_choice value items first
str3 = cell(1,numel(citems));
for k = 1:numel(citems)
    str3{k} = sprintf('* %s', subsref(citems{k}, substruct('.','name')));
end;
valitem = subsref(item, substruct('.','val'));
if ~isempty(valitem)
    str4 = cell(1,numel(valitem)+1);
    str4{1} = sprintf('Currently selected options:');
    for k = 1:numel(valitem)
        str4{k+1} = sprintf('* "%s"', subsref(valitem{k}, substruct('.','name')));
    end;
else
    str4 = {};
end;
str = [str1(:); {''}; str2(:); str3(:); str4(:)];