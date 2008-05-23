function dep = update_deps(dep, oid, nid)

% function dep = update_deps(dep, oid, nid)
% go through an array of dependencies and update tgt_exbranch and src_exbranch entries.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: update_deps.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

for k = 1:numel(dep)
    for l = 1:numel(oid)
        if isequal(dep(k).tgt_exbranch, oid{l})
            dep(k).tgt_exbranch = nid{l};
        end;
        if isequal(dep(k).src_exbranch, oid{l})
            dep(k).src_exbranch = nid{l};
        end;
    end;
end;