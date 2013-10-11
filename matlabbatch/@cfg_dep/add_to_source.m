function [cj, tdeps, cflag, dflag] = add_to_source(tdeps, cj)

% Add foreign target dependencies to own source dependencies
% function [cj tdeps cflag dflag] = add_to_source(tdeps, cj)
% This function adds target dependencies to the corresponding source
% exbranches. If a source exbranch can not be found at the exact location
% specified (e.g. because it has moved to a different level of the
% configuration hierarchy), it will try to find the corresponding exbranch
% and update the dependencies accordingly. If the dependencies can not be
% found, they will be marked for deletion. Note that update and deletion of
% dependencies for the current target item has to be done in the calling
% cfg_exbranch/harvest call.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: add_to_source.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

found = false(size(tdeps));
cflag = false(size(tdeps));
for k = 1:numel(tdeps)
    cflag(k) = false;
    try
        sitem    = subsref(cj, tdeps(k).src_exbranch); % Source item to deal with
        found(k) = isa(sitem, 'cfg_exbranch');
    catch
        found(k) = false;
    end;
    if ~found(k)
        try
            % look for correct source location
            citem = subsref(cj, tdeps(k).src_exbranch(1:2));
            id    = list(citem, cfg_findspec({{'class', 'cfg_exbranch'}}), cfg_tropts(cfg_findspec, 0, Inf, 0, 1, false));
            if ~isempty(id)
                otdep = tdeps(k);
                tdeps(k).src_exbranch = [tdeps(k).src_exbranch(1:2) id{1}];
                sitem    = subsref(cj, tdeps(k).src_exbranch); % Source item to deal with
                found(k) = isa(sitem, 'cfg_exbranch');
                cflag(k) = true;
            end;
        end;
    end;
    if found(k)
        if isempty(sitem.sdeps)
            sitem.sdeps = tdeps(k);
        else
            sitem.sdeps = [sitem.sdeps(:)' tdeps(k)];
        end;
        cj = subsasgn(cj, tdeps(k).src_exbranch,sitem);
        if cflag(k)
            % Update all references to a source once a change has been
            % detected
            cj = update_deps(cj, {otdep.src_exbranch}, {tdeps(k).src_exbranch});
        end
    end;
end;
tdeps = tdeps(found);
dflag = ~found;
cflag = cflag(found);
