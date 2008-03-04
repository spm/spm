function tn = tagnames(item, dflag)

% Return the tags of all children in the job tree of an item. dflag
% indicates whether the filled (false) or defaults (true) part of the
% tree should be searched. This function is identical for all in-tree
% nodes and not defined for leaf nodes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: tagnames.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

tp = treepart(item, dflag);
citems = subsref(item, substruct('.', tp));
tn = cell(size(citems));
for k = 1:numel(citems)
    tn{k} = gettag(citems{k});
end;