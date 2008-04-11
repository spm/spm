function tn = tagnames(item, dflag)

% function tn = tagnames(item, dflag)
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
% $Id: tagnames.m 1366 2008-04-11 10:24:17Z volkmar $

rev = '$Rev: 1366 $';

tp = treepart(item, dflag);
citems = subsref(item, substruct('.', tp));
tn = cell(size(citems));
for k = 1:numel(citems)
    tn{k} = gettag(citems{k});
end;