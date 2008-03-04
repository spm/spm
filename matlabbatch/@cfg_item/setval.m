function item = setval(item, val)

% function item = setval(item, val)
% set item.val{1} to val. Validity checks are performed through subsasgn.
% If val == {}, set item.val to {}. 
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: setval.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

if iscell(val) && isempty(val)
    item.val = {};
else
    item = subsasgn(item, substruct('.','val', '{}',{1}), val);
end;
