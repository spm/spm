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
% $Id: showmydoc.m 2085 2008-09-12 10:26:59Z volkmar $

rev = '$Rev: 2085 $'; %#ok
str1 = showdoc(item.cfg_item, indent);
citems = subsref(item, substruct('.','val'));
str2{1} = sprintf('This branch contains %d items:', numel(citems));
% Display short listing of branch items first
str3 = cell(1,numel(citems));
for k = 1:numel(citems)
    str3{k} = sprintf('* %s', subsref(citems{k}, substruct('.','name')));
end;
str = [str1(:); {''}; str2(:); str3(:)];