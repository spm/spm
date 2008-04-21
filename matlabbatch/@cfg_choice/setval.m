function item = setval(item, val)

% function item = setval(item, val)
% Set item.val{1} to item.values{val(1)}. If isempty(val), set item.val to {}.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: setval.m 1456 2008-04-21 15:03:41Z volkmar $

rev = '$Rev: 1456 $';

if isempty(val)
    item = subsasgn(item, substruct('.','val'), {});
else
    val = item.values{val(1)};
    item = subsasgn(item, substruct('.','val', '{}',{1}), val);
end;
