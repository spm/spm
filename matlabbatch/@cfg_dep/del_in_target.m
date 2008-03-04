function cj = del_in_target(sdeps, cj)

% If a dependency source has changed, drop all dependent (target)
% references recursively.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: del_in_target.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

% first, delete all immediate dependencies
for k = 1:numel(sdeps)
    eitem = subsref(cj, sdeps(k).tgt_exbranch); % dependent exbranch
    ditem = subsref(eitem, sdeps(k).tgt_input); % dependent item
    dind  = false(1,numel(ditem.val{1}));
    for l = 1:numel(ditem.val{1})
        dind(l) = ~isequalsource(sdeps(k),ditem.val{1}(l)); % ?.src_exbranch?
    end;
    if any(dind)
        ditem.val{1} = ditem.val{1}(dind); % Keep other dependencies
    else
        % nothing left
        ditem.val = {};
    end;
    eitem = subsasgn(eitem, sdeps(k).tgt_input, ditem);
    cj = subsasgn(cj, sdeps(k).tgt_exbranch, eitem);
end;
for k = 1:numel(sdeps)
    % re-harvest eitem - this may happen multiple times
    % better check for unique tgt_exbranch
    eitem = subsref(cj, sdeps(k).tgt_exbranch);
    [u1 u2 u3 u4 u5 cj] = harvest(eitem, cj, false, true);
end;