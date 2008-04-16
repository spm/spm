function cj = del_in_source(tdeps, cj)

% function cj = del_in_source(tdeps, cj)
% delete foreign target dependencies from own source dependencies.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: del_in_source.m 1431 2008-04-16 08:35:08Z volkmar $

rev = '$Rev: 1431 $';

for k = 1:numel(tdeps)
    sitem = subsref(cj, tdeps(k).src_exbranch); % Source item to deal with
    stind = false(1,numel(sitem.sdeps));
    for l = 1:numel(sitem.sdeps)
        stind(l) = ~isequaltarget(tdeps(k), sitem.sdeps(l));
    end;
    sitem.sdeps = sitem.sdeps(stind);
    cj = subsasgn(cj, tdeps(k).src_exbranch,sitem);
end;
