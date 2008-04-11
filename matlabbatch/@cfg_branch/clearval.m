function item = clearval(item, dflag)

% function item = clearval(item, dflag)
% Clear val fields in all items found in item.val.
% dflag is ignored in a cfg_branch.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: clearval.m 1373 2008-04-11 14:24:03Z spm $

rev = '$Rev: 1373 $';

for k = 1:numel(item.cfg_item.val)
    item.cfg_item.val{k} = clearval(item.cfg_item.val{k}, ...
                      dflag);
end;

